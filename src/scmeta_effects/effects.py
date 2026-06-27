from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy import sparse


class EffectError(RuntimeError):
    """Raised when effect-table inputs are inconsistent."""


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise EffectError(f"Expected JSON object: {path}")
    return payload


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _read_groups(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_pseudobulk(path: Path) -> tuple[sparse.csr_matrix, list[str], np.ndarray]:
    with h5py.File(path, "r") as handle:
        x = handle["X"]
        shape = tuple(int(value) for value in x.attrs["shape"])
        matrix = sparse.csr_matrix(
            (x["data"][()], x["indices"][()], x["indptr"][()]), shape=shape
        )
        feature_ids = [
            value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value)
            for value in handle["feature_id"][()]
        ]
        cell_counts = handle["cell_count"][()].astype(np.int64)
    return matrix, feature_ids, cell_counts


def _normalize_log1p(counts: np.ndarray, *, scale: float) -> np.ndarray:
    library_size = float(counts.sum())
    if library_size <= 0:
        return np.zeros_like(counts, dtype=np.float64)
    return np.log1p((counts / library_size) * scale)


def _safe_float(value: np.ndarray, index: int) -> float:
    return float(value[index])


def _write_effect_h5(
    *,
    path: Path,
    dataset_id: str,
    feature_ids: list[str],
    contrast_rows: list[dict[str, Any]],
    delta: np.ndarray,
    treated_log: np.ndarray,
    control_log: np.ndarray,
) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    if temporary.exists():
        temporary.unlink()
    with h5py.File(temporary, "w") as handle:
        handle.attrs["dataset_id"] = dataset_id
        handle.attrs["encoding-type"] = "effect_matrix"
        handle.attrs["normalization"] = "log1p(CPM)"
        handle.create_dataset("delta_log1p_cpm", data=delta.astype(np.float32), compression="gzip", shuffle=True)
        handle.create_dataset("treated_log1p_cpm", data=treated_log.astype(np.float32), compression="gzip", shuffle=True)
        handle.create_dataset("control_log1p_cpm", data=control_log.astype(np.float32), compression="gzip", shuffle=True)
        handle.create_dataset(
            "feature_id",
            data=np.asarray(feature_ids, dtype=h5py.string_dtype("utf-8")),
            compression="gzip",
        )
        handle.create_dataset(
            "contrast_id",
            data=np.asarray([row["contrast_id"] for row in contrast_rows], dtype=h5py.string_dtype("utf-8")),
            compression="gzip",
        )
    temporary.replace(path)
    return {
        "path": str(path),
        "contrasts": len(contrast_rows),
        "features": len(feature_ids),
        "sha256": _sha256(path),
    }


def _write_contrasts(rows: list[dict[str, Any]], path: Path) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "contrast_id",
        "dataset_id",
        "treated_group_id",
        "control_group_ids",
        "perturbation_id",
        "perturbation_type",
        "perturbation_target",
        "guide_id",
        "drug_id",
        "treated_cell_count",
        "control_cell_count",
        "treated_library_size",
        "control_library_size",
        "status",
        "reason",
        "cell_line",
        "cell_type",
        "disease",
        "dose",
        "dose_unit",
        "time",
        "time_unit",
        "treated_condition_ids",
        "control_condition_ids",
        "treated_condition_count",
        "control_condition_count",
        "matched_strata_count",
        "model",
    ]
    temporary = path.with_suffix(path.suffix + ".partial")
    with path.with_suffix(path.suffix + ".partial").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    temporary.replace(path)
    return {"path": str(path), "rows": len(rows), "sha256": _sha256(path)}


CONDITION_FIELDS = (
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


def _condition_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(field, "") for field in CONDITION_FIELDS)


def _split_ids(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def _write_top_genes(
    *,
    path: Path,
    dataset_id: str,
    feature_ids: list[str],
    contrast_rows: list[dict[str, Any]],
    delta: np.ndarray,
    treated_log: np.ndarray,
    control_log: np.ndarray,
    top_n: int,
) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    rows_written = 0
    with gzip.open(temporary, "wt", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "dataset_id",
            "contrast_id",
            "perturbation_id",
            "gene_id",
            "effect_size",
            "treated_log1p_cpm",
            "control_log1p_cpm",
            "rank_abs_effect",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for contrast_index, contrast in enumerate(contrast_rows):
            values = delta[contrast_index]
            n = min(top_n, len(feature_ids))
            if n == 0:
                continue
            selected = np.argpartition(np.abs(values), -n)[-n:]
            selected = selected[np.argsort(-np.abs(values[selected]))]
            for rank, gene_index in enumerate(selected, start=1):
                writer.writerow(
                    {
                        "dataset_id": dataset_id,
                        "contrast_id": contrast["contrast_id"],
                        "perturbation_id": contrast["perturbation_id"],
                        "gene_id": feature_ids[int(gene_index)],
                        "effect_size": _safe_float(values, int(gene_index)),
                        "treated_log1p_cpm": _safe_float(treated_log[contrast_index], int(gene_index)),
                        "control_log1p_cpm": _safe_float(control_log[contrast_index], int(gene_index)),
                        "rank_abs_effect": rank,
                    }
                )
                rows_written += 1
    temporary.replace(path)
    return {"path": str(path), "rows": rows_written, "sha256": _sha256(path)}


def _dataset_effects(
    *,
    dataset: dict[str, Any],
    output_dir: Path,
    scale: float,
    min_cells: int,
    top_genes: int,
) -> dict[str, Any]:
    dataset_id = dataset["dataset_id"]
    groups = _read_groups(_resolve_path(dataset["groups"]["path"]))
    matrix, feature_ids, cell_counts = _read_pseudobulk(_resolve_path(dataset["matrix"]["path"]))
    if len(groups) != matrix.shape[0]:
        raise EffectError(f"Group count does not match matrix rows for {dataset_id}")

    control_indices = [i for i, row in enumerate(groups) if row.get("is_control") == "true"]
    treated_indices = [
        i
        for i, row in enumerate(groups)
        if row.get("is_control") != "true" and int(row.get("cell_count", 0)) >= min_cells
    ]
    if not control_indices:
        raise EffectError(f"No control groups found for {dataset_id}")
    if not treated_indices:
        raise EffectError(f"No treated groups pass min cell threshold for {dataset_id}")

    control_counts = np.asarray(matrix[control_indices].sum(axis=0)).ravel()
    control_cell_count = int(cell_counts[control_indices].sum())
    control_library_size = float(control_counts.sum())
    control_log = _normalize_log1p(control_counts, scale=scale)

    contrast_rows = []
    treated_logs = []
    control_logs = []
    deltas = []
    for contrast_number, group_index in enumerate(treated_indices, start=1):
        treated_counts = np.asarray(matrix[group_index].todense()).ravel()
        treated_log = _normalize_log1p(treated_counts, scale=scale)
        row = groups[group_index]
        contrast_rows.append(
            {
                "contrast_id": f"{dataset_id}_contrast_{contrast_number:07d}",
                "dataset_id": dataset_id,
                "treated_group_id": row["group_id"],
                "control_group_ids": ";".join(groups[index]["group_id"] for index in control_indices),
                "perturbation_id": row.get("perturbation_id", ""),
                "perturbation_type": row.get("perturbation_type", ""),
                "perturbation_target": row.get("perturbation_target", ""),
                "guide_id": row.get("guide_id", ""),
                "drug_id": row.get("drug_id", ""),
                "treated_cell_count": int(cell_counts[group_index]),
                "control_cell_count": control_cell_count,
                "treated_library_size": float(treated_counts.sum()),
                "control_library_size": control_library_size,
                "model": "delta_log1p_cpm_intervention_pseudobulk",
            }
        )
        treated_logs.append(treated_log)
        control_logs.append(control_log)
        deltas.append(treated_log - control_log)

    treated_matrix = np.vstack(treated_logs)
    control_matrix = np.vstack(control_logs)
    delta_matrix = np.vstack(deltas)
    dataset_dir = output_dir / dataset_id
    effect_matrix = _write_effect_h5(
        path=dataset_dir / "effects.h5",
        dataset_id=dataset_id,
        feature_ids=feature_ids,
        contrast_rows=contrast_rows,
        delta=delta_matrix,
        treated_log=treated_matrix,
        control_log=control_matrix,
    )
    contrasts = _write_contrasts(contrast_rows, dataset_dir / "contrasts.csv")
    top = _write_top_genes(
        path=dataset_dir / "top_genes.csv.gz",
        dataset_id=dataset_id,
        feature_ids=feature_ids,
        contrast_rows=contrast_rows,
        delta=delta_matrix,
        treated_log=treated_matrix,
        control_log=control_matrix,
        top_n=top_genes,
    )
    return {
        "dataset_id": dataset_id,
        "controls": len(control_indices),
        "treated_groups_input": len([row for row in groups if row.get("is_control") != "true"]),
        "contrasts": len(contrast_rows),
        "min_cells": min_cells,
        "matrix": effect_matrix,
        "contrast_table": contrasts,
        "top_genes": top,
    }


def _dataset_pooled_effects(
    *,
    dataset: dict[str, Any],
    conditions: dict[str, dict[str, str]],
    pooled_contrasts: list[dict[str, str]],
    output_dir: Path,
    scale: float,
    top_genes: int,
) -> dict[str, Any]:
    dataset_id = dataset["dataset_id"]
    groups = _read_groups(_resolve_path(dataset["groups"]["path"]))
    matrix, feature_ids, cell_counts = _read_pseudobulk(_resolve_path(dataset["matrix"]["path"]))
    if len(groups) != matrix.shape[0]:
        raise EffectError(f"Group count does not match matrix rows for {dataset_id}")

    group_by = dataset.get("group_by", "")
    if group_by != "condition":
        raise EffectError(
            f"Pooled effects require condition-level pseudobulk for {dataset_id}; got {group_by!r}"
        )

    group_index_by_condition = {_condition_key(row): index for index, row in enumerate(groups)}
    condition_index_by_id = {}
    for condition_id, row in conditions.items():
        if row.get("dataset_id") != dataset_id:
            continue
        index = group_index_by_condition.get(_condition_key(row))
        if index is not None:
            condition_index_by_id[condition_id] = index

    contrast_rows = []
    treated_logs = []
    control_logs = []
    deltas = []
    skipped = 0
    for contrast in pooled_contrasts:
        if contrast.get("dataset_id") != dataset_id:
            continue
        treated_indices = [
            condition_index_by_id[condition_id]
            for condition_id in _split_ids(contrast.get("treated_condition_ids", ""))
            if condition_id in condition_index_by_id
        ]
        control_indices = [
            condition_index_by_id[condition_id]
            for condition_id in _split_ids(contrast.get("control_condition_ids", ""))
            if condition_id in condition_index_by_id
        ]
        if not treated_indices or not control_indices:
            skipped += 1
            continue

        treated_counts = np.asarray(matrix[treated_indices].sum(axis=0)).ravel()
        control_counts = np.asarray(matrix[control_indices].sum(axis=0)).ravel()
        treated_log = _normalize_log1p(treated_counts, scale=scale)
        control_log = _normalize_log1p(control_counts, scale=scale)
        treated_library_size = float(treated_counts.sum())
        control_library_size = float(control_counts.sum())
        contrast_rows.append(
            {
                "contrast_id": contrast["contrast_id"],
                "dataset_id": dataset_id,
                "treated_group_id": ";".join(groups[index]["group_id"] for index in treated_indices),
                "control_group_ids": ";".join(groups[index]["group_id"] for index in control_indices),
                "perturbation_id": contrast.get("perturbation_id", ""),
                "treated_cell_count": int(cell_counts[treated_indices].sum()),
                "control_cell_count": int(cell_counts[control_indices].sum()),
                "treated_library_size": treated_library_size,
                "control_library_size": control_library_size,
                "status": contrast.get("status", ""),
                "reason": contrast.get("reason", ""),
                "cell_line": contrast.get("cell_line", ""),
                "cell_type": contrast.get("cell_type", ""),
                "disease": contrast.get("disease", ""),
                "dose": contrast.get("dose", ""),
                "dose_unit": contrast.get("dose_unit", ""),
                "time": contrast.get("time", ""),
                "time_unit": contrast.get("time_unit", ""),
                "treated_condition_ids": contrast.get("treated_condition_ids", ""),
                "control_condition_ids": contrast.get("control_condition_ids", ""),
                "treated_condition_count": contrast.get("treated_condition_count", ""),
                "control_condition_count": contrast.get("control_condition_count", ""),
                "matched_strata_count": contrast.get("matched_strata_count", ""),
                "model": "delta_log1p_cpm_pooled_condition_pseudobulk",
            }
        )
        treated_logs.append(treated_log)
        control_logs.append(control_log)
        deltas.append(treated_log - control_log)

    if contrast_rows:
        treated_matrix = np.vstack(treated_logs)
        control_matrix = np.vstack(control_logs)
        delta_matrix = np.vstack(deltas)
    else:
        treated_matrix = np.empty((0, len(feature_ids)), dtype=np.float64)
        control_matrix = np.empty((0, len(feature_ids)), dtype=np.float64)
        delta_matrix = np.empty((0, len(feature_ids)), dtype=np.float64)

    dataset_dir = output_dir / dataset_id
    effect_matrix = _write_effect_h5(
        path=dataset_dir / "effects.h5",
        dataset_id=dataset_id,
        feature_ids=feature_ids,
        contrast_rows=contrast_rows,
        delta=delta_matrix,
        treated_log=treated_matrix,
        control_log=control_matrix,
    )
    contrasts = _write_contrasts(contrast_rows, dataset_dir / "contrasts.csv")
    top = _write_top_genes(
        path=dataset_dir / "top_genes.csv.gz",
        dataset_id=dataset_id,
        feature_ids=feature_ids,
        contrast_rows=contrast_rows,
        delta=delta_matrix,
        treated_log=treated_matrix,
        control_log=control_matrix,
        top_n=top_genes,
    )
    return {
        "dataset_id": dataset_id,
        "contrasts": len(contrast_rows),
        "skipped_contrasts": skipped,
        "matrix": effect_matrix,
        "contrast_table": contrasts,
        "top_genes": top,
    }


def compute_effects(
    *,
    pseudobulk_manifest_path: Path,
    output_dir: Path,
    include_datasets: set[str] | None = None,
    scale: float = 1_000_000.0,
    min_cells: int = 50,
    top_genes: int = 100,
) -> dict[str, Any]:
    manifest = _read_json(pseudobulk_manifest_path)
    if "datasets" not in manifest or not isinstance(manifest["datasets"], list):
        raise EffectError("Pseudobulk manifest must contain a datasets list")
    outputs = []
    for dataset in manifest["datasets"]:
        dataset_id = str(dataset["dataset_id"])
        if include_datasets and dataset_id not in include_datasets:
            continue
        outputs.append(
            _dataset_effects(
                dataset=dataset,
                output_dir=output_dir,
                scale=scale,
                min_cells=min_cells,
                top_genes=top_genes,
            )
        )
    result = {
        "schema_version": "1.0",
        "release_id": manifest.get("release_id"),
        "pseudobulk_manifest": str(pseudobulk_manifest_path),
        "normalization": "log1p(CPM)",
        "model": "treated_minus_control_delta",
        "scale": scale,
        "min_cells": min_cells,
        "top_genes": top_genes,
        "datasets": outputs,
    }
    _write_json(result, output_dir / "effects-manifest.json")
    return result


def compute_pooled_effects(
    *,
    pseudobulk_manifest_path: Path,
    pooled_contrasts_path: Path,
    output_dir: Path,
    include_statuses: set[str] | None = None,
    include_datasets: set[str] | None = None,
    scale: float = 1_000_000.0,
    top_genes: int = 100,
) -> dict[str, Any]:
    statuses = include_statuses or {"pass", "warn"}
    manifest = _read_json(pseudobulk_manifest_path)
    if manifest.get("group_by") != "condition":
        raise EffectError("Pooled effects require a condition-level pseudobulk manifest")
    harmonization_manifest_path = _resolve_path(manifest["harmonization_manifest"])
    harmonization = _read_json(harmonization_manifest_path)
    conditions_path = _resolve_path(harmonization["conditions"]["path"])
    conditions = {
        row["condition_id"]: row
        for row in _read_csv(conditions_path)
        if row.get("condition_id")
    }
    pooled_contrasts = [
        row
        for row in _read_csv(pooled_contrasts_path)
        if row.get("status") in statuses
    ]

    outputs = []
    for dataset in manifest["datasets"]:
        dataset_id = str(dataset["dataset_id"])
        if include_datasets and dataset_id not in include_datasets:
            continue
        outputs.append(
            _dataset_pooled_effects(
                dataset=dataset,
                conditions=conditions,
                pooled_contrasts=pooled_contrasts,
                output_dir=output_dir,
                scale=scale,
                top_genes=top_genes,
            )
        )
    result = {
        "schema_version": "1.0",
        "release_id": manifest.get("release_id"),
        "pseudobulk_manifest": str(pseudobulk_manifest_path),
        "pooled_contrasts": str(pooled_contrasts_path),
        "normalization": "log1p(CPM)",
        "model": "treated_minus_matched_control_delta",
        "scale": scale,
        "top_genes": top_genes,
        "included_statuses": sorted(statuses),
        "input_contrasts": len(pooled_contrasts),
        "datasets": outputs,
    }
    _write_json(result, output_dir / "effects-manifest.json")
    return result
