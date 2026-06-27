from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from scmeta_curate.h5ad import decode_array


class QCError(RuntimeError):
    """Raised when QC inputs are missing or inconsistent."""


PERTURBATION_COLUMNS = ("perturbation", "condition")
CONTROL_STATUS_COLUMNS = ("control_status", "is_control")
BATCH_COLUMNS = ("batch", "batch_id", "replicate", "replicate_id", "plate", "gem_group")
CONTROL_VALUES = {
    "control",
    "negative_control",
    "non-targeting",
    "nontargeting",
    "non_targeting",
    "scramble",
    "vehicle",
    "dmso",
    "untreated",
    "mock",
}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise QCError(f"Expected JSON object: {path}")
    return payload


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _write_csv(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if row.get(key) is None else row.get(key) for key in fieldnames})
    temporary.replace(path)


def _resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _profile_for_dataset(record: dict[str, Any]) -> dict[str, Any]:
    profile_path = record.get("lineage", {}).get("profile")
    if not profile_path:
        raise QCError(f"Release dataset has no profile path: {record.get('dataset_id')}")
    path = _resolve_path(profile_path)
    if not path.exists():
        raise QCError(f"Profile file does not exist: {path}")
    return _read_json(path)


def _first_existing(group: h5py.Group, names: tuple[str, ...]) -> str | None:
    for name in names:
        if name in group:
            return name
    return None


def _missing_value(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and np.isnan(value):
        return True
    text = str(value).strip()
    return text == "" or text.lower() in {"nan", "none", "null", "na", "<missing>"}


def _normalize(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip().lower().replace("_", "-")


def _control_mask(
    perturbations: np.ndarray,
    *,
    control_status: np.ndarray | None = None,
) -> np.ndarray:
    mask = np.zeros(len(perturbations), dtype=bool)
    if control_status is not None:
        for index, value in enumerate(control_status):
            normalized = _normalize(value)
            if normalized in CONTROL_VALUES or normalized in {"true", "1"}:
                mask[index] = True
    for index, value in enumerate(perturbations):
        normalized = _normalize(value)
        if normalized in CONTROL_VALUES:
            mask[index] = True
    return mask


def _status_from_checks(failures: list[str], warnings: list[str]) -> str:
    if failures:
        return "fail"
    if warnings:
        return "warn"
    return "pass"


def _deduplicate(values: list[str]) -> list[str]:
    seen = set()
    result = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def _load_obs_vectors(path: Path) -> dict[str, np.ndarray | str | None]:
    with h5py.File(path, "r") as h5ad:
        if "obs" not in h5ad:
            raise QCError(f"H5AD is missing obs: {path}")
        obs = h5ad["obs"]
        perturbation_column = _first_existing(obs, PERTURBATION_COLUMNS)
        if perturbation_column is None:
            raise QCError(f"H5AD is missing perturbation column: {path}")
        control_column = _first_existing(obs, CONTROL_STATUS_COLUMNS)
        batch_column = _first_existing(obs, BATCH_COLUMNS)
        return {
            "perturbation_column": perturbation_column,
            "control_column": control_column,
            "batch_column": batch_column,
            "perturbation": decode_array(obs, perturbation_column),
            "control_status": None
            if control_column is None
            else decode_array(obs, control_column),
            "batch": None if batch_column is None else decode_array(obs, batch_column),
        }


def _dataset_qc(
    record: dict[str, Any],
    profile: dict[str, Any],
    obs: dict[str, np.ndarray | str | None],
    *,
    min_control_cells: int,
) -> dict[str, Any]:
    dataset_id = record["dataset_id"]
    path = _resolve_path(record["analysis_input"]["path"])
    perturbations = obs["perturbation"]
    control_status = obs["control_status"]
    assert isinstance(perturbations, np.ndarray)
    assert control_status is None or isinstance(control_status, np.ndarray)
    controls = _control_mask(perturbations, control_status=control_status)
    missing_perturbations = int(sum(_missing_value(value) for value in perturbations))

    failures: list[str] = []
    warnings: list[str] = []
    if not path.exists():
        failures.append("analysis input file is missing")
    if record.get("admission_status", "").startswith("pass") is False:
        failures.append(f"release admission is not pass-like: {record.get('admission_status')}")
    if not profile.get("matrix", {}).get("integer_like_sample"):
        failures.append("sampled matrix values are not integer-like")
    if profile.get("matrix", {}).get("n_obs") != len(perturbations):
        failures.append("profile n_obs does not match observed perturbation vector")
    if missing_perturbations:
        failures.append(f"missing perturbation labels: {missing_perturbations}")
    if int(controls.sum()) < min_control_cells:
        failures.append(f"control cells below threshold: {int(controls.sum())} < {min_control_cells}")

    if not profile.get("counts", {}).get("raw_present") and not profile.get("counts", {}).get("layers"):
        warnings.append("no raw object or count layer; X is treated as counts")
    if obs["batch_column"] is None:
        warnings.append("no explicit batch/replicate column found")
    warnings.extend(profile.get("warnings", []))
    warnings = _deduplicate(warnings)

    return {
        "dataset_id": dataset_id,
        "status": _status_from_checks(failures, warnings),
        "admission_status": record.get("admission_status"),
        "analysis_input": str(path),
        "n_obs": len(perturbations),
        "n_vars": profile.get("matrix", {}).get("n_vars"),
        "matrix_encoding": profile.get("matrix", {}).get("encoding"),
        "matrix_integer_like_sample": profile.get("matrix", {}).get("integer_like_sample"),
        "counts_location": profile.get("counts", {}).get("counts_location_inferred"),
        "perturbation_column": obs["perturbation_column"],
        "control_column": obs["control_column"],
        "batch_column": obs["batch_column"],
        "perturbation_count": len(set(str(value) for value in perturbations if not _missing_value(value))),
        "missing_perturbation_cells": missing_perturbations,
        "control_cells": int(controls.sum()),
        "failures": failures,
        "warnings": warnings,
    }


def _perturbation_qc(
    dataset_id: str,
    obs: dict[str, np.ndarray | str | None],
    *,
    min_cells: int,
    min_control_cells: int,
) -> list[dict[str, Any]]:
    perturbations = obs["perturbation"]
    control_status = obs["control_status"]
    batches = obs["batch"]
    assert isinstance(perturbations, np.ndarray)
    assert control_status is None or isinstance(control_status, np.ndarray)
    assert batches is None or isinstance(batches, np.ndarray)
    controls = _control_mask(perturbations, control_status=control_status)
    total_controls = int(controls.sum())

    counts: Counter[str] = Counter()
    control_counts: Counter[str] = Counter()
    batch_sets: dict[str, set[str]] = defaultdict(set)
    for index, value in enumerate(perturbations):
        perturbation = "<missing>" if _missing_value(value) else str(value)
        counts[perturbation] += 1
        if controls[index]:
            control_counts[perturbation] += 1
        if batches is not None and not _missing_value(batches[index]):
            batch_sets[perturbation].add(str(batches[index]))

    rows = []
    for perturbation, cell_count in sorted(counts.items()):
        is_control = control_counts[perturbation] == cell_count
        failures = []
        warnings = []
        if perturbation == "<missing>":
            failures.append("missing perturbation label")
        if not is_control and cell_count < min_cells:
            failures.append(f"cell count below threshold: {cell_count} < {min_cells}")
        if not is_control and total_controls < min_control_cells:
            failures.append(f"dataset control cells below threshold: {total_controls} < {min_control_cells}")
        if batches is None:
            warnings.append("batch/replicate confounding cannot be assessed")
        elif not is_control and len(batch_sets[perturbation]) <= 1:
            warnings.append("perturbation observed in one batch/replicate level")
        rows.append(
            {
                "dataset_id": dataset_id,
                "perturbation_id": perturbation,
                "status": "control" if is_control else _status_from_checks(failures, warnings),
                "is_control": is_control,
                "cell_count": cell_count,
                "control_cells_in_dataset": total_controls,
                "batch_column": obs["batch_column"],
                "batch_level_count": None if batches is None else len(batch_sets[perturbation]),
                "failures": "; ".join(failures),
                "warnings": "; ".join(warnings),
            }
        )
    return rows


def run_qc(
    *,
    release_manifest_path: Path,
    output_dir: Path,
    min_cells: int = 50,
    min_control_cells: int = 50,
) -> dict[str, Any]:
    release = _read_json(release_manifest_path)
    if "datasets" not in release or not isinstance(release["datasets"], list):
        raise QCError("Release manifest must contain a datasets list")

    dataset_rows = []
    perturbation_rows = []
    for record in release["datasets"]:
        profile = _profile_for_dataset(record)
        h5ad_path = _resolve_path(record["analysis_input"]["path"])
        obs = _load_obs_vectors(h5ad_path)
        dataset_rows.append(
            _dataset_qc(record, profile, obs, min_control_cells=min_control_cells)
        )
        perturbation_rows.extend(
            _perturbation_qc(
                record["dataset_id"],
                obs,
                min_cells=min_cells,
                min_control_cells=min_control_cells,
            )
        )

    dataset_counts = Counter(row["status"] for row in dataset_rows)
    perturbation_counts = Counter(row["status"] for row in perturbation_rows)
    summary = {
        "schema_version": "1.0",
        "release_id": release.get("release_id"),
        "release_manifest": str(release_manifest_path),
        "thresholds": {
            "min_cells_per_perturbation": min_cells,
            "min_control_cells_per_dataset": min_control_cells,
        },
        "control_detection": {
            "control_status_columns": list(CONTROL_STATUS_COLUMNS),
            "perturbation_columns": list(PERTURBATION_COLUMNS),
            "control_values": sorted(CONTROL_VALUES),
        },
        "dataset_count": len(dataset_rows),
        "dataset_status_counts": dict(dataset_counts),
        "perturbation_count": len(perturbation_rows),
        "perturbation_status_counts": dict(perturbation_counts),
        "outputs": {
            "dataset_qc_json": str(output_dir / "dataset-qc.json"),
            "dataset_qc_csv": str(output_dir / "dataset-qc.csv"),
            "perturbation_qc_json": str(output_dir / "perturbation-qc.json"),
            "perturbation_qc_csv": str(output_dir / "perturbation-qc.csv"),
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_json({"schema_version": "1.0", "datasets": dataset_rows}, output_dir / "dataset-qc.json")
    _write_csv(dataset_rows, output_dir / "dataset-qc.csv")
    _write_json(
        {"schema_version": "1.0", "perturbations": perturbation_rows},
        output_dir / "perturbation-qc.json",
    )
    _write_csv(perturbation_rows, output_dir / "perturbation-qc.csv")
    _write_json(summary, output_dir / "qc-summary.json")
    return summary
