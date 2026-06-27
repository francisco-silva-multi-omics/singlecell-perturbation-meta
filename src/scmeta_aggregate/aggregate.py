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
from scipy import sparse


class AggregateError(RuntimeError):
    """Raised when pseudobulk aggregation inputs are inconsistent."""


GROUP_FIELDS = {
    "intervention": (
        "dataset_id",
        "perturbation_id",
        "perturbation_type",
        "perturbation_target",
        "guide_id",
        "drug_id",
        "is_control",
    ),
    "bio-condition": (
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
    ),
    "condition": (
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
    ),
}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise AggregateError(f"Expected JSON object: {path}")
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


def _read_gzip_csv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _read_feature_ids(path: Path) -> list[str]:
    rows = _read_gzip_csv(path)
    ids = []
    for row in rows:
        ids.append(row.get("feature_id") or row.get("gene_id") or row.get("gene_symbol") or "")
    return ids


def _group_cells(
    obs_path: Path,
    *,
    dataset_id: str,
    group_by: str,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    fields = GROUP_FIELDS[group_by]
    group_lookup: dict[tuple[str, ...], int] = {}
    group_ids = []
    group_rows: list[dict[str, Any]] = []
    counts: Counter[int] = Counter()

    with gzip.open(obs_path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = tuple(row.get(field, "") for field in fields)
            group_index = group_lookup.get(key)
            if group_index is None:
                group_index = len(group_rows)
                group_lookup[key] = group_index
                group_id = f"{dataset_id}_{group_by.replace('-', '_')}_{group_index + 1:07d}"
                metadata = {"group_id": group_id, "group_by": group_by}
                metadata.update({field: value for field, value in zip(fields, key)})
                group_rows.append(metadata)
            group_ids.append(group_index)
            counts[group_index] += 1

    for index, row in enumerate(group_rows):
        row["cell_count"] = counts[index]
    return np.asarray(group_ids, dtype=np.int32), group_rows


def _write_group_metadata(rows: list[dict[str, Any]], path: Path) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row.keys()})
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    temporary.replace(path)
    return {"path": str(path), "rows": len(rows), "sha256": _sha256(path)}


def _write_sparse_h5(
    matrix: sparse.spmatrix,
    *,
    output_path: Path,
    dataset_id: str,
    group_by: str,
    group_rows: list[dict[str, Any]],
    feature_ids: list[str],
) -> dict[str, Any]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".partial")
    if temporary.exists():
        temporary.unlink()
    csr = matrix.tocsr()
    with h5py.File(temporary, "w") as handle:
        handle.attrs["dataset_id"] = dataset_id
        handle.attrs["group_by"] = group_by
        handle.attrs["encoding-type"] = "pseudobulk_csr_matrix"
        handle.attrs["shape"] = np.asarray(csr.shape, dtype=np.int64)
        x = handle.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["shape"] = np.asarray(csr.shape, dtype=np.int64)
        x.create_dataset("data", data=csr.data, compression="gzip", shuffle=True)
        x.create_dataset("indices", data=csr.indices, compression="gzip", shuffle=True)
        x.create_dataset("indptr", data=csr.indptr, compression="gzip", shuffle=True)
        handle.create_dataset(
            "group_id",
            data=np.asarray([row["group_id"] for row in group_rows], dtype=h5py.string_dtype("utf-8")),
            compression="gzip",
        )
        handle.create_dataset(
            "feature_id",
            data=np.asarray(feature_ids, dtype=h5py.string_dtype("utf-8")),
            compression="gzip",
        )
        handle.create_dataset(
            "cell_count",
            data=np.asarray([row["cell_count"] for row in group_rows], dtype=np.int64),
            compression="gzip",
        )
    temporary.replace(output_path)
    return {
        "path": str(output_path),
        "rows": int(csr.shape[0]),
        "columns": int(csr.shape[1]),
        "nnz": int(csr.nnz),
        "sha256": _sha256(output_path),
    }


def _csr_chunk_matrix(x: h5py.Group, start: int, end: int, n_vars: int) -> sparse.csr_matrix:
    indptr = x["indptr"][start : end + 1]
    data_start = int(indptr[0])
    data_end = int(indptr[-1])
    chunk_indptr = indptr - data_start
    data = x["data"][data_start:data_end]
    indices = x["indices"][data_start:data_end]
    return sparse.csr_matrix((data, indices, chunk_indptr), shape=(end - start, n_vars))


def _group_indicator(group_ids: np.ndarray, n_groups: int) -> sparse.csr_matrix:
    cols = np.arange(len(group_ids), dtype=np.int32)
    data = np.ones(len(group_ids), dtype=np.float64)
    return sparse.csr_matrix((data, (group_ids, cols)), shape=(n_groups, len(group_ids)))


def _aggregate_csr(
    x: h5py.Group,
    *,
    group_ids: np.ndarray,
    n_groups: int,
    chunk_size: int,
) -> sparse.csr_matrix:
    shape = tuple(int(value) for value in x.attrs["shape"])
    n_obs, n_vars = shape
    if n_obs != len(group_ids):
        raise AggregateError("group id count does not match X rows")
    result = sparse.csr_matrix((n_groups, n_vars), dtype=np.float64)
    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        chunk = _csr_chunk_matrix(x, start, end, n_vars)
        indicator = _group_indicator(group_ids[start:end], n_groups)
        result = result + indicator @ chunk
    return result


def _aggregate_csc(
    x: h5py.Group,
    *,
    group_ids: np.ndarray,
    n_groups: int,
) -> sparse.csc_matrix:
    shape = tuple(int(value) for value in x.attrs["shape"])
    n_obs, n_vars = shape
    if n_obs != len(group_ids):
        raise AggregateError("group id count does not match X rows")
    indptr = x["indptr"][()]
    out_indptr = np.empty(n_vars + 1, dtype=np.int64)
    out_indptr[0] = 0
    out_indices = []
    out_data = []
    nnz = 0
    for column in range(n_vars):
        start, end = int(indptr[column]), int(indptr[column + 1])
        rows = x["indices"][start:end]
        values = x["data"][start:end]
        if len(rows):
            sums = np.bincount(group_ids[rows], weights=values, minlength=n_groups)
            groups = np.flatnonzero(sums)
            out_indices.append(groups.astype(np.int32))
            out_data.append(sums[groups])
            nnz += len(groups)
        out_indptr[column + 1] = nnz
    indices = np.concatenate(out_indices) if out_indices else np.asarray([], dtype=np.int32)
    data = np.concatenate(out_data) if out_data else np.asarray([], dtype=np.float64)
    return sparse.csc_matrix((data, indices, out_indptr), shape=(n_groups, n_vars))


def _aggregate_dense(
    x: h5py.Dataset,
    *,
    group_ids: np.ndarray,
    n_groups: int,
    chunk_size: int,
) -> sparse.csr_matrix:
    n_obs, n_vars = x.shape
    if n_obs != len(group_ids):
        raise AggregateError("group id count does not match X rows")
    result = np.zeros((n_groups, n_vars), dtype=np.float64)
    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        chunk = np.asarray(x[start:end], dtype=np.float64)
        indicator = _group_indicator(group_ids[start:end], n_groups)
        result += indicator @ chunk
    return sparse.csr_matrix(result)


def _aggregate_matrix(
    h5ad_path: Path,
    *,
    group_ids: np.ndarray,
    n_groups: int,
    chunk_size: int,
) -> sparse.spmatrix:
    with h5py.File(h5ad_path, "r") as handle:
        x = handle["X"]
        if isinstance(x, h5py.Dataset):
            return _aggregate_dense(
                x,
                group_ids=group_ids,
                n_groups=n_groups,
                chunk_size=chunk_size,
            )
        encoding = x.attrs.get("encoding-type")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding == "csr_matrix":
            return _aggregate_csr(
                x,
                group_ids=group_ids,
                n_groups=n_groups,
                chunk_size=chunk_size,
            )
        if encoding == "csc_matrix":
            return _aggregate_csc(
                x,
                group_ids=group_ids,
                n_groups=n_groups,
            )
        raise AggregateError(f"Unsupported X encoding: {encoding}")


def aggregate_harmonized(
    *,
    harmonization_manifest_path: Path,
    output_dir: Path,
    group_by: str = "intervention",
    include_datasets: set[str] | None = None,
    chunk_size: int = 4096,
) -> dict[str, Any]:
    if group_by not in GROUP_FIELDS:
        raise AggregateError(f"Unsupported group-by value: {group_by}")
    manifest = _read_json(harmonization_manifest_path)
    if "datasets" not in manifest or not isinstance(manifest["datasets"], list):
        raise AggregateError("Harmonization manifest must contain a datasets list")

    outputs = []
    for dataset in manifest["datasets"]:
        dataset_id = str(dataset["dataset_id"])
        if include_datasets and dataset_id not in include_datasets:
            continue
        h5ad_path = _resolve_path(dataset["analysis_input"])
        obs_path = _resolve_path(dataset["obs"]["path"])
        var_path = _resolve_path(dataset["var"]["path"])
        group_ids, group_rows = _group_cells(obs_path, dataset_id=dataset_id, group_by=group_by)
        feature_ids = _read_feature_ids(var_path)
        matrix = _aggregate_matrix(
            h5ad_path,
            group_ids=group_ids,
            n_groups=len(group_rows),
            chunk_size=chunk_size,
        )
        dataset_dir = output_dir / group_by / dataset_id
        group_metadata = _write_group_metadata(group_rows, dataset_dir / "groups.csv")
        matrix_output = _write_sparse_h5(
            matrix,
            output_path=dataset_dir / "pseudobulk.h5",
            dataset_id=dataset_id,
            group_by=group_by,
            group_rows=group_rows,
            feature_ids=feature_ids,
        )
        outputs.append(
            {
                "dataset_id": dataset_id,
                "analysis_input": str(h5ad_path),
                "group_by": group_by,
                "groups": group_metadata,
                "matrix": matrix_output,
            }
        )

    aggregate_manifest = {
        "schema_version": "1.0",
        "release_id": manifest.get("release_id"),
        "harmonization_manifest": str(harmonization_manifest_path),
        "group_by": group_by,
        "chunk_size": chunk_size,
        "datasets": outputs,
    }
    _write_json(aggregate_manifest, output_dir / group_by / "pseudobulk-manifest.json")
    return aggregate_manifest
