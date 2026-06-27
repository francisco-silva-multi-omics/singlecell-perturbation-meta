from __future__ import annotations

from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from scmeta_curate.h5ad import dataframe_index, decode_array


class HarmonizeError(RuntimeError):
    """Raised when harmonization inputs are inconsistent."""


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

OBS_FIELDS = [
    "dataset_id",
    "cell_id",
    "cell_line",
    "cell_type",
    "disease",
    "organism",
    "perturbation_id",
    "perturbation_type",
    "perturbation_target",
    "guide_id",
    "drug_id",
    "dose",
    "dose_unit",
    "time",
    "time_unit",
    "control_status",
    "is_control",
    "batch_id",
    "replicate_id",
    "n_counts",
    "n_genes",
    "percent_mito",
    "percent_ribo",
    "source_cell_barcode",
]

VAR_FIELDS = [
    "dataset_id",
    "feature_id",
    "gene_id",
    "gene_symbol",
    "ensembl_id",
    "chromosome",
    "start",
    "end",
    "gene_type",
    "highly_variable",
]

OBS_ALIASES = {
    "cell_line": ("cell_line",),
    "cell_type": ("cell_type", "celltype"),
    "disease": ("disease", "disease_status"),
    "organism": ("organism", "species"),
    "perturbation_id": ("perturbation", "condition"),
    "perturbation_type": ("perturbation_type",),
    "perturbation_target": ("perturbation_target", "target", "gene"),
    "guide_id": ("guide_id",),
    "drug_id": ("drug_id", "drug", "drug_name"),
    "dose": ("dose", "dose_value", "drug_dose"),
    "dose_unit": ("dose_unit", "drug_dose_unit"),
    "time": ("time", "drug_time", "timepoint"),
    "time_unit": ("time_unit", "drug_time_unit"),
    "control_status": ("control_status",),
    "batch_id": ("batch_id", "batch", "plate"),
    "replicate_id": ("replicate_id", "replicate"),
    "n_counts": ("n_counts", "ncounts", "UMI_count", "UMI count"),
    "n_genes": ("n_genes", "ngenes"),
    "percent_mito": ("percent_mito",),
    "percent_ribo": ("percent_ribo",),
    "source_cell_barcode": ("source_cell_barcode", "cell_barcode"),
}

VAR_ALIASES = {
    "gene_id": ("gene_id",),
    "gene_symbol": ("gene_symbol", "gene_name"),
    "ensembl_id": ("ensembl_id",),
    "chromosome": ("chromosome", "chr"),
    "start": ("start",),
    "end": ("end",),
    "gene_type": ("gene_type", "class"),
    "highly_variable": ("highly_variable",),
}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise HarmonizeError(f"Expected JSON object: {path}")
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _stringify(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and np.isnan(value):
        return ""
    return str(value)


def _normalize(value: Any) -> str:
    return _stringify(value).strip().lower().replace("_", "-")


def _is_control(perturbation: Any, control_status: Any = None) -> bool:
    status = _normalize(control_status)
    if status in CONTROL_VALUES or status in {"true", "1"}:
        return True
    return _normalize(perturbation) in CONTROL_VALUES


def _control_status(perturbation: Any, control_status: Any = None) -> str:
    if _stringify(control_status):
        return _stringify(control_status)
    return "control" if _is_control(perturbation, control_status) else "targeting"


def _first_existing(group: h5py.Group, names: tuple[str, ...]) -> str | None:
    for name in names:
        if name in group:
            return name
    return None


def _read_optional_column(group: h5py.Group, aliases: tuple[str, ...]) -> tuple[str | None, np.ndarray | None]:
    name = _first_existing(group, aliases)
    if name is None:
        return None, None
    return name, decode_array(group, name)


def _constant(length: int, value: str = "") -> np.ndarray:
    return np.full(length, value, dtype=object)


def _load_dataframe_columns(
    group: h5py.Group,
    aliases: dict[str, tuple[str, ...]],
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, str | None]]:
    index = dataframe_index(group)
    length = len(index)
    values = {}
    sources = {}
    for field, field_aliases in aliases.items():
        source, column = _read_optional_column(group, field_aliases)
        sources[field] = source
        values[field] = _constant(length) if column is None else column
    return index, values, sources


def _open_gzip_csv(path: Path, fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    handle = gzip.open(temporary, "wt", newline="", encoding="utf-8")
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    return temporary, handle, writer


def _promote(temporary: Path, output: Path) -> None:
    temporary.replace(output)


def _write_plain_csv(rows: list[dict[str, Any]], path: Path, fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    temporary.replace(path)


def _obs_row(
    *,
    dataset_id: str,
    cell_id: Any,
    columns: dict[str, np.ndarray],
    index: int,
) -> dict[str, str]:
    perturbation = columns["perturbation_id"][index]
    control_status = columns["control_status"][index]
    perturbation_type = _stringify(columns["perturbation_type"][index])
    drug_id = _stringify(columns["drug_id"][index])
    if not drug_id and perturbation_type.lower() == "drug":
        drug_id = _stringify(perturbation)
    source_cell_barcode = _stringify(columns["source_cell_barcode"][index]) or _stringify(cell_id)
    return {
        "dataset_id": dataset_id,
        "cell_id": _stringify(cell_id),
        "cell_line": _stringify(columns["cell_line"][index]),
        "cell_type": _stringify(columns["cell_type"][index]),
        "disease": _stringify(columns["disease"][index]),
        "organism": _stringify(columns["organism"][index]),
        "perturbation_id": _stringify(perturbation),
        "perturbation_type": perturbation_type,
        "perturbation_target": _stringify(columns["perturbation_target"][index]),
        "guide_id": _stringify(columns["guide_id"][index]),
        "drug_id": drug_id,
        "dose": _stringify(columns["dose"][index]),
        "dose_unit": _stringify(columns["dose_unit"][index]),
        "time": _stringify(columns["time"][index]),
        "time_unit": _stringify(columns["time_unit"][index]),
        "control_status": _control_status(perturbation, control_status),
        "is_control": str(_is_control(perturbation, control_status)).lower(),
        "batch_id": _stringify(columns["batch_id"][index]),
        "replicate_id": _stringify(columns["replicate_id"][index]),
        "n_counts": _stringify(columns["n_counts"][index]),
        "n_genes": _stringify(columns["n_genes"][index]),
        "percent_mito": _stringify(columns["percent_mito"][index]),
        "percent_ribo": _stringify(columns["percent_ribo"][index]),
        "source_cell_barcode": source_cell_barcode,
    }


def _condition_key(row: dict[str, str]) -> tuple[str, ...]:
    fields = (
        "dataset_id",
        "cell_line",
        "cell_type",
        "disease",
        "perturbation_id",
        "dose",
        "dose_unit",
        "time",
        "time_unit",
        "control_status",
        "is_control",
        "batch_id",
        "replicate_id",
    )
    return tuple(row[field] for field in fields)


def _intervention_key(row: dict[str, str]) -> tuple[str, ...]:
    fields = (
        "dataset_id",
        "perturbation_id",
        "perturbation_type",
        "perturbation_target",
        "guide_id",
        "drug_id",
        "is_control",
    )
    return tuple(row[field] for field in fields)


def _write_obs(
    h5ad: h5py.File,
    *,
    dataset_id: str,
    output_path: Path,
) -> tuple[dict[str, Any], Counter[tuple[str, ...]], Counter[tuple[str, ...]]]:
    index, columns, sources = _load_dataframe_columns(h5ad["obs"], OBS_ALIASES)
    temporary, handle, writer = _open_gzip_csv(output_path, OBS_FIELDS)
    conditions: Counter[tuple[str, ...]] = Counter()
    interventions: Counter[tuple[str, ...]] = Counter()
    try:
        for row_index, cell_id in enumerate(index):
            row = _obs_row(
                dataset_id=dataset_id,
                cell_id=cell_id,
                columns=columns,
                index=row_index,
            )
            writer.writerow(row)
            conditions[_condition_key(row)] += 1
            interventions[_intervention_key(row)] += 1
    finally:
        handle.close()
    _promote(temporary, output_path)
    return {
        "path": str(output_path),
        "rows": int(len(index)),
        "sha256": _sha256(output_path),
        "column_sources": sources,
    }, conditions, interventions


def _write_var(
    h5ad: h5py.File,
    *,
    dataset_id: str,
    output_path: Path,
) -> dict[str, Any]:
    index, columns, sources = _load_dataframe_columns(h5ad["var"], VAR_ALIASES)
    temporary, handle, writer = _open_gzip_csv(output_path, VAR_FIELDS)
    try:
        for row_index, feature_id in enumerate(index):
            row = {
                "dataset_id": dataset_id,
                "feature_id": _stringify(feature_id),
                "gene_id": _stringify(columns["gene_id"][row_index]) or _stringify(feature_id),
                "gene_symbol": _stringify(columns["gene_symbol"][row_index]) or _stringify(feature_id),
                "ensembl_id": _stringify(columns["ensembl_id"][row_index]),
                "chromosome": _stringify(columns["chromosome"][row_index]),
                "start": _stringify(columns["start"][row_index]),
                "end": _stringify(columns["end"][row_index]),
                "gene_type": _stringify(columns["gene_type"][row_index]),
                "highly_variable": _stringify(columns["highly_variable"][row_index]),
            }
            writer.writerow(row)
    finally:
        handle.close()
    _promote(temporary, output_path)
    return {
        "path": str(output_path),
        "rows": int(len(index)),
        "sha256": _sha256(output_path),
        "column_sources": sources,
    }


def _condition_rows(counter: Counter[tuple[str, ...]]) -> list[dict[str, Any]]:
    fields = (
        "dataset_id",
        "cell_line",
        "cell_type",
        "disease",
        "perturbation_id",
        "dose",
        "dose_unit",
        "time",
        "time_unit",
        "control_status",
        "is_control",
        "batch_id",
        "replicate_id",
    )
    rows = []
    for number, (key, count) in enumerate(sorted(counter.items()), start=1):
        row = dict(zip(fields, key))
        row["condition_id"] = f"condition_{number:07d}"
        row["cell_count"] = count
        rows.append(row)
    return rows


def _intervention_rows(counter: Counter[tuple[str, ...]]) -> list[dict[str, Any]]:
    fields = (
        "dataset_id",
        "perturbation_id",
        "perturbation_type",
        "perturbation_target",
        "guide_id",
        "drug_id",
        "is_control",
    )
    rows = []
    for number, (key, count) in enumerate(sorted(counter.items()), start=1):
        row = dict(zip(fields, key))
        row["intervention_id"] = f"intervention_{number:07d}"
        row["cell_count"] = count
        rows.append(row)
    return rows


def harmonize_release(
    *,
    release_manifest_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    release = _read_json(release_manifest_path)
    if "datasets" not in release or not isinstance(release["datasets"], list):
        raise HarmonizeError("Release manifest must contain a datasets list")

    obs_dir = output_dir / "obs"
    var_dir = output_dir / "var"
    dataset_outputs = []
    all_conditions: Counter[tuple[str, ...]] = Counter()
    all_interventions: Counter[tuple[str, ...]] = Counter()

    for record in release["datasets"]:
        dataset_id = str(record["dataset_id"])
        h5ad_path = _resolve_path(record["analysis_input"]["path"])
        if not h5ad_path.exists():
            raise HarmonizeError(f"Analysis input does not exist: {h5ad_path}")
        obs_path = obs_dir / f"{dataset_id}.csv.gz"
        var_path = var_dir / f"{dataset_id}.csv.gz"
        with h5py.File(h5ad_path, "r") as h5ad:
            if "obs" not in h5ad or "var" not in h5ad:
                raise HarmonizeError(f"H5AD must contain obs and var: {h5ad_path}")
            obs_output, conditions, interventions = _write_obs(
                h5ad,
                dataset_id=dataset_id,
                output_path=obs_path,
            )
            var_output = _write_var(
                h5ad,
                dataset_id=dataset_id,
                output_path=var_path,
            )
        all_conditions.update(conditions)
        all_interventions.update(interventions)
        dataset_outputs.append(
            {
                "dataset_id": dataset_id,
                "analysis_input": str(h5ad_path),
                "obs": obs_output,
                "var": var_output,
            }
        )

    condition_rows = _condition_rows(all_conditions)
    intervention_rows = _intervention_rows(all_interventions)
    conditions_path = output_dir / "conditions.csv"
    interventions_path = output_dir / "interventions.csv"
    _write_plain_csv(
        intervention_rows,
        interventions_path,
        [
            "intervention_id",
            "dataset_id",
            "perturbation_id",
            "perturbation_type",
            "perturbation_target",
            "guide_id",
            "drug_id",
            "is_control",
            "cell_count",
        ],
    )
    _write_plain_csv(
        condition_rows,
        conditions_path,
        [
            "condition_id",
            "dataset_id",
            "cell_line",
            "cell_type",
            "disease",
            "perturbation_id",
            "dose",
            "dose_unit",
            "time",
            "time_unit",
            "control_status",
            "is_control",
            "batch_id",
            "replicate_id",
            "cell_count",
        ],
    )

    manifest = {
        "schema_version": "1.0",
        "release_id": release.get("release_id"),
        "release_manifest": str(release_manifest_path),
        "output_dir": str(output_dir),
        "datasets": dataset_outputs,
        "obs_schema": OBS_FIELDS,
        "var_schema": VAR_FIELDS,
        "interventions": {
            "path": str(interventions_path),
            "rows": len(intervention_rows),
            "sha256": _sha256(interventions_path),
        },
        "conditions": {
            "path": str(conditions_path),
            "rows": len(condition_rows),
            "sha256": _sha256(conditions_path),
        },
    }
    _write_json(manifest, output_dir / "harmonization-manifest.json")
    return manifest
