from __future__ import annotations

from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from typing import Any


class ContrastError(RuntimeError):
    """Raised when contrast inputs are missing or inconsistent."""


STRICT_FIELDS = (
    "dataset_id",
    "cell_line",
    "cell_type",
    "disease",
    "time",
    "time_unit",
    "batch_id",
    "replicate_id",
)

BIO_FIELDS = (
    "dataset_id",
    "cell_line",
    "cell_type",
    "disease",
    "time",
    "time_unit",
)

DATASET_FIELDS = ("dataset_id",)

POOLED_CONTEXT_FIELDS = (
    "dataset_id",
    "cell_line",
    "cell_type",
    "disease",
    "perturbation_id",
    "dose",
    "dose_unit",
    "time",
    "time_unit",
)

POOLED_CONTROL_FIELDS = (
    "dataset_id",
    "cell_line",
    "cell_type",
    "disease",
    "time",
    "time_unit",
    "batch_id",
    "replicate_id",
)

POOLED_STRATUM_FIELDS = ("batch_id", "replicate_id")


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ContrastError(f"Expected JSON object: {path}")
    return payload


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _read_conditions(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(rows: list[dict[str, Any]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    temporary.replace(path)


def _key(row: dict[str, str], fields: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(row.get(field, "") for field in fields)


def _index_controls(
    controls: list[dict[str, str]], fields: tuple[str, ...]
) -> dict[tuple[str, ...], list[dict[str, str]]]:
    index: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in controls:
        index[_key(row, fields)].append(row)
    return index


def _group_rows(
    rows: list[dict[str, str]], fields: tuple[str, ...]
) -> dict[tuple[str, ...], list[dict[str, str]]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[_key(row, fields)].append(row)
    return grouped


def _match_controls(
    treated: dict[str, str],
    *,
    strict: dict[tuple[str, ...], list[dict[str, str]]],
    bio: dict[tuple[str, ...], list[dict[str, str]]],
    dataset: dict[tuple[str, ...], list[dict[str, str]]],
) -> tuple[str, list[dict[str, str]]]:
    matches = strict.get(_key(treated, STRICT_FIELDS), [])
    if matches:
        return "strict_context_batch_replicate", matches
    matches = bio.get(_key(treated, BIO_FIELDS), [])
    if matches:
        return "bio_context", matches
    matches = dataset.get(_key(treated, DATASET_FIELDS), [])
    if matches:
        return "dataset_fallback", matches
    return "unmatched", []


def _admission_status(
    *,
    match_rule: str,
    treated_cells: int,
    control_cells: int,
    min_treated_cells: int,
    min_control_cells: int,
) -> str:
    if match_rule == "unmatched":
        return "fail"
    if treated_cells < min_treated_cells or control_cells < min_control_cells:
        return "fail"
    if match_rule == "dataset_fallback":
        return "warn"
    return "pass"


def _reason(
    *,
    match_rule: str,
    treated_cells: int,
    control_cells: int,
    min_treated_cells: int,
    min_control_cells: int,
) -> str:
    reasons = []
    if match_rule == "unmatched":
        reasons.append("no matched controls")
    if treated_cells < min_treated_cells:
        reasons.append(f"treated cells below threshold: {treated_cells} < {min_treated_cells}")
    if control_cells < min_control_cells:
        reasons.append(f"control cells below threshold: {control_cells} < {min_control_cells}")
    if match_rule == "dataset_fallback":
        reasons.append("controls matched only at dataset level")
    return "; ".join(reasons)


def build_contrasts(
    *,
    harmonization_manifest_path: Path,
    output_dir: Path,
    min_treated_cells: int = 50,
    min_control_cells: int = 50,
) -> dict[str, Any]:
    manifest = _read_json(harmonization_manifest_path)
    conditions_path = _resolve_path(manifest["conditions"]["path"])
    conditions = _read_conditions(conditions_path)
    controls = [row for row in conditions if row.get("is_control") == "true"]
    treated = [row for row in conditions if row.get("is_control") != "true"]
    strict = _index_controls(controls, STRICT_FIELDS)
    bio = _index_controls(controls, BIO_FIELDS)
    dataset = _index_controls(controls, DATASET_FIELDS)

    rows = []
    for number, treated_row in enumerate(treated, start=1):
        match_rule, matches = _match_controls(
            treated_row,
            strict=strict,
            bio=bio,
            dataset=dataset,
        )
        treated_cells = int(treated_row.get("cell_count") or 0)
        control_cells = sum(int(row.get("cell_count") or 0) for row in matches)
        status = _admission_status(
            match_rule=match_rule,
            treated_cells=treated_cells,
            control_cells=control_cells,
            min_treated_cells=min_treated_cells,
            min_control_cells=min_control_cells,
        )
        rows.append(
            {
                "contrast_id": f"condition_contrast_{number:08d}",
                "dataset_id": treated_row.get("dataset_id", ""),
                "treated_condition_id": treated_row.get("condition_id", ""),
                "control_condition_ids": ";".join(row.get("condition_id", "") for row in matches),
                "match_rule": match_rule,
                "status": status,
                "reason": _reason(
                    match_rule=match_rule,
                    treated_cells=treated_cells,
                    control_cells=control_cells,
                    min_treated_cells=min_treated_cells,
                    min_control_cells=min_control_cells,
                ),
                "perturbation_id": treated_row.get("perturbation_id", ""),
                "cell_line": treated_row.get("cell_line", ""),
                "cell_type": treated_row.get("cell_type", ""),
                "disease": treated_row.get("disease", ""),
                "dose": treated_row.get("dose", ""),
                "dose_unit": treated_row.get("dose_unit", ""),
                "time": treated_row.get("time", ""),
                "time_unit": treated_row.get("time_unit", ""),
                "batch_id": treated_row.get("batch_id", ""),
                "replicate_id": treated_row.get("replicate_id", ""),
                "treated_cell_count": treated_cells,
                "control_cell_count": control_cells,
                "control_match_count": len(matches),
            }
        )

    fieldnames = [
        "contrast_id",
        "dataset_id",
        "treated_condition_id",
        "control_condition_ids",
        "match_rule",
        "status",
        "reason",
        "perturbation_id",
        "cell_line",
        "cell_type",
        "disease",
        "dose",
        "dose_unit",
        "time",
        "time_unit",
        "batch_id",
        "replicate_id",
        "treated_cell_count",
        "control_cell_count",
        "control_match_count",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    contrast_path = output_dir / "condition-contrasts.csv"
    _write_csv(rows, contrast_path, fieldnames)

    status_counts = Counter(row["status"] for row in rows)
    rule_counts = Counter(row["match_rule"] for row in rows)
    dataset_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        dataset_counts[row["dataset_id"]][row["status"]] += 1
    summary = {
        "schema_version": "1.0",
        "release_id": manifest.get("release_id"),
        "harmonization_manifest": str(harmonization_manifest_path),
        "conditions": str(conditions_path),
        "output": str(contrast_path),
        "thresholds": {
            "min_treated_cells": min_treated_cells,
            "min_control_cells": min_control_cells,
        },
        "matching_rules": {
            "strict_context_batch_replicate": list(STRICT_FIELDS),
            "bio_context": list(BIO_FIELDS),
            "dataset_fallback": list(DATASET_FIELDS),
        },
        "treated_condition_count": len(treated),
        "control_condition_count": len(controls),
        "contrast_count": len(rows),
        "status_counts": dict(status_counts),
        "match_rule_counts": dict(rule_counts),
        "dataset_status_counts": {
            dataset_id: dict(counts) for dataset_id, counts in sorted(dataset_counts.items())
        },
    }
    _write_json(summary, output_dir / "contrast-summary.json")
    return summary


def _pooled_status(
    *,
    treated_cells: int,
    control_cells: int,
    matched_strata: int,
    min_treated_cells: int,
    min_control_cells: int,
    min_strata: int,
) -> str:
    if matched_strata == 0:
        return "fail"
    if treated_cells < min_treated_cells or control_cells < min_control_cells:
        return "fail"
    if matched_strata < min_strata:
        return "warn"
    return "pass"


def _pooled_reason(
    *,
    treated_cells: int,
    control_cells: int,
    matched_strata: int,
    min_treated_cells: int,
    min_control_cells: int,
    min_strata: int,
) -> str:
    reasons = []
    if matched_strata == 0:
        reasons.append("no matched control strata")
    if treated_cells < min_treated_cells:
        reasons.append(f"treated cells below pooled threshold: {treated_cells} < {min_treated_cells}")
    if control_cells < min_control_cells:
        reasons.append(f"control cells below pooled threshold: {control_cells} < {min_control_cells}")
    if 0 < matched_strata < min_strata:
        reasons.append(f"matched strata below recommendation: {matched_strata} < {min_strata}")
    return "; ".join(reasons)


def build_pooled_contrasts(
    *,
    harmonization_manifest_path: Path,
    output_dir: Path,
    min_treated_cells: int = 50,
    min_control_cells: int = 50,
    min_strata: int = 2,
) -> dict[str, Any]:
    manifest = _read_json(harmonization_manifest_path)
    conditions_path = _resolve_path(manifest["conditions"]["path"])
    conditions = _read_conditions(conditions_path)
    controls = [row for row in conditions if row.get("is_control") == "true"]
    treated = [row for row in conditions if row.get("is_control") != "true"]
    control_by_stratum = _index_controls(controls, POOLED_CONTROL_FIELDS)
    treated_groups = _group_rows(treated, POOLED_CONTEXT_FIELDS)

    rows = []
    for number, (context_key, treated_rows) in enumerate(sorted(treated_groups.items()), start=1):
        context = dict(zip(POOLED_CONTEXT_FIELDS, context_key))
        matched_controls = []
        matched_treated = []
        for row in treated_rows:
            control_key = _key(row, POOLED_CONTROL_FIELDS)
            controls_for_stratum = control_by_stratum.get(control_key, [])
            if controls_for_stratum:
                matched_treated.append(row)
                matched_controls.extend(controls_for_stratum)
        unique_control_ids = {
            row.get("condition_id", ""): row for row in matched_controls
        }
        strata = {
            _key(row, POOLED_STRATUM_FIELDS)
            for row in matched_treated
        }
        treated_cells = sum(int(row.get("cell_count") or 0) for row in matched_treated)
        control_cells = sum(
            int(row.get("cell_count") or 0) for row in unique_control_ids.values()
        )
        matched_strata = len(strata)
        status = _pooled_status(
            treated_cells=treated_cells,
            control_cells=control_cells,
            matched_strata=matched_strata,
            min_treated_cells=min_treated_cells,
            min_control_cells=min_control_cells,
            min_strata=min_strata,
        )
        rows.append(
            {
                "contrast_id": f"pooled_contrast_{number:08d}",
                "dataset_id": context["dataset_id"],
                "perturbation_id": context["perturbation_id"],
                "dose": context["dose"],
                "dose_unit": context["dose_unit"],
                "time": context["time"],
                "time_unit": context["time_unit"],
                "cell_line": context["cell_line"],
                "cell_type": context["cell_type"],
                "disease": context["disease"],
                "status": status,
                "reason": _pooled_reason(
                    treated_cells=treated_cells,
                    control_cells=control_cells,
                    matched_strata=matched_strata,
                    min_treated_cells=min_treated_cells,
                    min_control_cells=min_control_cells,
                    min_strata=min_strata,
                ),
                "match_rule": "pooled_context_strata",
                "treated_condition_ids": ";".join(row.get("condition_id", "") for row in matched_treated),
                "control_condition_ids": ";".join(sorted(unique_control_ids)),
                "treated_cell_count": treated_cells,
                "control_cell_count": control_cells,
                "treated_condition_count": len(matched_treated),
                "control_condition_count": len(unique_control_ids),
                "matched_strata_count": matched_strata,
            }
        )

    fieldnames = [
        "contrast_id",
        "dataset_id",
        "perturbation_id",
        "dose",
        "dose_unit",
        "time",
        "time_unit",
        "cell_line",
        "cell_type",
        "disease",
        "status",
        "reason",
        "match_rule",
        "treated_condition_ids",
        "control_condition_ids",
        "treated_cell_count",
        "control_cell_count",
        "treated_condition_count",
        "control_condition_count",
        "matched_strata_count",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    contrast_path = output_dir / "pooled-contrasts.csv"
    _write_csv(rows, contrast_path, fieldnames)

    status_counts = Counter(row["status"] for row in rows)
    dataset_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        dataset_counts[row["dataset_id"]][row["status"]] += 1
    summary = {
        "schema_version": "1.0",
        "release_id": manifest.get("release_id"),
        "harmonization_manifest": str(harmonization_manifest_path),
        "conditions": str(conditions_path),
        "output": str(contrast_path),
        "thresholds": {
            "min_treated_cells": min_treated_cells,
            "min_control_cells": min_control_cells,
            "min_strata": min_strata,
        },
        "matching_rule": {
            "pooled_context": list(POOLED_CONTEXT_FIELDS),
            "control_strata": list(POOLED_CONTROL_FIELDS),
            "stratum_identity": list(POOLED_STRATUM_FIELDS),
        },
        "contrast_count": len(rows),
        "status_counts": dict(status_counts),
        "dataset_status_counts": {
            dataset_id: dict(counts) for dataset_id, counts in sorted(dataset_counts.items())
        },
    }
    _write_json(summary, output_dir / "pooled-contrast-summary.json")
    return summary
