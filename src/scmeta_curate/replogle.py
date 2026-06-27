from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from scmeta_curate.h5ad import CurationError, dataframe_index, decode_array


def _logical_matrix_digest(dataset: h5py.Dataset, rows_per_chunk: int = 256) -> str:
    digest = hashlib.sha256()
    for start in range(0, dataset.shape[0], rows_per_chunk):
        block = np.ascontiguousarray(dataset[start : start + rows_per_chunk])
        digest.update(block.tobytes())
    return digest.hexdigest()


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def verify_replogle(primary_path: Path, harmonized_path: Path) -> dict:
    with h5py.File(primary_path, "r") as primary, h5py.File(harmonized_path, "r") as harmonized:
        if not isinstance(primary["X"], h5py.Dataset) or not isinstance(
            harmonized["X"], h5py.Dataset
        ):
            raise CurationError("Replogle verification expects dense X datasets")
        shape_match = primary["X"].shape == harmonized["X"].shape
        matrix_shape = list(primary["X"].shape)
        if not shape_match:
            raise CurationError(
                f"Matrix shape mismatch: {primary['X'].shape} != {harmonized['X'].shape}"
            )

        primary_cells = dataframe_index(primary["obs"])
        harmonized_cells = dataframe_index(harmonized["obs"])
        cell_match = np.array_equal(primary_cells, harmonized_cells)

        primary_gene_ids = dataframe_index(primary["var"])
        harmonized_gene_ids = decode_array(harmonized["var"], "ensembl_id")
        gene_match = np.array_equal(primary_gene_ids, harmonized_gene_ids)

        primary_digest = _logical_matrix_digest(primary["X"])
        harmonized_digest = _logical_matrix_digest(harmonized["X"])
        matrix_match = primary_digest == harmonized_digest

    return {
        "status": "pass" if shape_match and cell_match and gene_match and matrix_match else "fail",
        "primary_path": str(primary_path),
        "harmonized_path": str(harmonized_path),
        "shape": matrix_shape,
        "shape_match": shape_match,
        "cell_ids_match": cell_match,
        "gene_ids_match": gene_match,
        "matrix_values_match": matrix_match,
        "primary_file_size_bytes": primary_path.stat().st_size,
        "primary_file_sha256": _file_digest(primary_path),
        "primary_matrix_sha256": primary_digest,
        "harmonized_matrix_sha256": harmonized_digest,
    }


def write_verification(result: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(output_path.suffix + ".partial")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output_path)
