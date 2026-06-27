from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_curate.h5ad import decode_array, filter_h5ad_rows
from scmeta_curate.replogle import verify_replogle


def _string_dataset(parent, name, values):
    dataset = parent.create_dataset(
        name, data=np.asarray(values, dtype=h5py.string_dtype("utf-8"))
    )
    dataset.attrs["encoding-type"] = "string-array"
    dataset.attrs["encoding-version"] = "0.2.0"


def _categorical(parent, name, values):
    categories = sorted(set(values))
    lookup = {value: index for index, value in enumerate(categories)}
    group = parent.create_group(name)
    group.attrs["encoding-type"] = "categorical"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["ordered"] = False
    _string_dataset(group, "categories", categories)
    codes = group.create_dataset(
        "codes", data=np.asarray([lookup[value] for value in values], dtype=np.int8)
    )
    codes.attrs["encoding-type"] = "array"
    codes.attrs["encoding-version"] = "0.2.0"


def _write_fixture(path: Path, matrix: np.ndarray, encoding: str = "csr_matrix"):
    with h5py.File(path, "w") as h5ad:
        h5ad.attrs["encoding-type"] = "anndata"
        h5ad.attrs["encoding-version"] = "0.1.0"
        x = h5ad.create_group("X")
        x.attrs["encoding-type"] = encoding
        x.attrs["encoding-version"] = "0.1.0"
        x.attrs["shape"] = matrix.shape
        if encoding == "csr_matrix":
            rows, columns = np.nonzero(matrix)
            order = np.lexsort((columns, rows))
            rows, columns = rows[order], columns[order]
            counts = np.bincount(rows, minlength=matrix.shape[0])
            indices = columns
        else:
            rows, columns = np.nonzero(matrix)
            order = np.lexsort((rows, columns))
            rows, columns = rows[order], columns[order]
            counts = np.bincount(columns, minlength=matrix.shape[1])
            indices = rows
        data = matrix[rows, columns]
        x.create_dataset("data", data=data)
        x.create_dataset("indices", data=indices.astype(np.int32))
        x.create_dataset(
            "indptr", data=np.concatenate([[0], np.cumsum(counts)]).astype(np.int32)
        )

        obs = h5ad.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "cell_id"
        obs.attrs["column-order"] = np.asarray(
            ["perturbation"], dtype=h5py.string_dtype("utf-8")
        )
        _string_dataset(obs, "cell_id", [f"cell-{i}" for i in range(matrix.shape[0])])
        _categorical(obs, "perturbation", ["a", "control", "b", "control"])

        var = h5ad.create_group("var")
        var.attrs["encoding-type"] = "dataframe"
        var.attrs["encoding-version"] = "0.2.0"
        var.attrs["_index"] = "gene_id"
        var.attrs["column-order"] = np.asarray([], dtype=h5py.string_dtype("utf-8"))
        _string_dataset(var, "gene_id", [f"gene-{i}" for i in range(matrix.shape[1])])
        for key in ("layers", "obsm", "obsp", "uns", "varm", "varp"):
            h5ad.create_group(key)


def _dense(group):
    shape = tuple(group.attrs["shape"])
    result = np.zeros(shape, dtype=group["data"].dtype)
    indptr = group["indptr"][()]
    if group.attrs["encoding-type"] == "csr_matrix":
        for row in range(shape[0]):
            first, last = indptr[row : row + 2]
            result[row, group["indices"][first:last]] = group["data"][first:last]
    else:
        for column in range(shape[1]):
            first, last = indptr[column : column + 2]
            result[group["indices"][first:last], column] = group["data"][first:last]
    return result


class CurationTests(unittest.TestCase):
    def test_filters_csr_rows_and_updates_obs(self):
        self._assert_sparse_filter("csr_matrix")

    def test_filters_csc_rows_and_updates_obs(self):
        self._assert_sparse_filter("csc_matrix")

    def _assert_sparse_filter(self, encoding):
        matrix = np.asarray([[1, 0, 2], [0, 3, 0], [4, 5, 0], [0, 0, 6]])
        keep = np.asarray([True, False, True, True])
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.h5ad"
            output = Path(temporary) / "output.h5ad"
            _write_fixture(source, matrix, encoding)
            filter_h5ad_rows(
                source,
                output,
                keep,
                obs_updates={"nperts": np.asarray([1, 0, 1, 0])},
            )
            with h5py.File(output, "r") as h5ad:
                np.testing.assert_array_equal(_dense(h5ad["X"]), matrix[keep])
                self.assertEqual(decode_array(h5ad["obs"], "cell_id").tolist(), ["cell-0", "cell-2", "cell-3"])
                self.assertEqual(decode_array(h5ad["obs"], "nperts").tolist(), [1, 1, 0])

    def test_verifies_equal_dense_replogle_matrices(self):
        matrix = np.arange(12, dtype=np.float32).reshape(4, 3)
        with tempfile.TemporaryDirectory() as temporary:
            primary = Path(temporary) / "primary.h5ad"
            harmonized = Path(temporary) / "harmonized.h5ad"
            for path, local in ((primary, False), (harmonized, True)):
                with h5py.File(path, "w") as h5ad:
                    h5ad.create_dataset("X", data=matrix)
                    obs = h5ad.create_group("obs")
                    obs.attrs["_index"] = "cell_id"
                    _string_dataset(obs, "cell_id", [f"c{i}" for i in range(4)])
                    var = h5ad.create_group("var")
                    var.attrs["_index"] = "feature_id"
                    _string_dataset(var, "feature_id", [f"g{i}" for i in range(3)])
                    if local:
                        _string_dataset(var, "ensembl_id", [f"g{i}" for i in range(3)])
            result = verify_replogle(primary, harmonized)
            self.assertEqual(result["status"], "pass")
            self.assertTrue(result["matrix_values_match"])


if __name__ == "__main__":
    unittest.main()
