from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_profile.profiler import (
    profile_from_curation_manifest,
    profile_from_release_manifest,
    profile_h5ad,
)


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


def _write_dense_fixture(path: Path):
    with h5py.File(path, "w") as h5ad:
        h5ad.attrs["encoding-type"] = "anndata"
        h5ad.create_dataset("X", data=np.asarray([[1, 0], [2, 3], [0, 4]]))

        obs = h5ad.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "cell_id"
        obs.attrs["column-order"] = np.asarray(
            ["perturbation", "batch"], dtype=h5py.string_dtype("utf-8")
        )
        _string_dataset(obs, "cell_id", ["c1", "c2", "c3"])
        _categorical(obs, "perturbation", ["control", "target", "target"])
        _categorical(obs, "batch", ["b1", "b1", "b2"])

        var = h5ad.create_group("var")
        var.attrs["encoding-type"] = "dataframe"
        var.attrs["encoding-version"] = "0.2.0"
        var.attrs["_index"] = "gene_id"
        var.attrs["column-order"] = np.asarray(
            ["gene_symbol"], dtype=h5py.string_dtype("utf-8")
        )
        _string_dataset(var, "gene_id", ["ENSG1", "ENSG2"])
        _string_dataset(var, "gene_symbol", ["A", "B"])
        for key in ("layers", "obsm", "obsp", "uns", "varm", "varp"):
            h5ad.create_group(key)


def _write_sparse_fixture(path: Path):
    with h5py.File(path, "w") as h5ad:
        x = h5ad.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["shape"] = np.asarray([2, 3])
        x.create_dataset("data", data=np.asarray([1, 2, 3], dtype=np.int32))
        x.create_dataset("indices", data=np.asarray([0, 2, 1], dtype=np.int32))
        x.create_dataset("indptr", data=np.asarray([0, 2, 3], dtype=np.int32))

        obs = h5ad.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        obs.attrs["column-order"] = np.asarray(["perturbation"], dtype=h5py.string_dtype("utf-8"))
        _string_dataset(obs, "cell_id", ["c1", "c2"])
        _categorical(obs, "perturbation", ["control", "gene"])

        var = h5ad.create_group("var")
        var.attrs["_index"] = "gene_id"
        var.attrs["column-order"] = np.asarray([], dtype=h5py.string_dtype("utf-8"))
        _string_dataset(var, "gene_id", ["g1", "g2", "g3"])


class ProfileTests(unittest.TestCase):
    def test_profiles_dense_h5ad(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "dense.h5ad"
            _write_dense_fixture(path)
            profile = profile_h5ad(path, dataset_id="fixture")
            self.assertEqual(profile["matrix"]["shape"], [3, 2])
            self.assertEqual(profile["matrix"]["encoding"], "dense")
            self.assertTrue(profile["matrix"]["integer_like_sample"])
            self.assertEqual(
                profile["obs"]["summaries"]["perturbation"]["top_values"][0],
                {"value": "target", "count": 2},
            )
            self.assertIn("no raw object or count layer found", "; ".join(profile["warnings"]))

    def test_profiles_sparse_h5ad(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "sparse.h5ad"
            _write_sparse_fixture(path)
            profile = profile_h5ad(path, dataset_id="sparse")
            self.assertEqual(profile["matrix"]["encoding"], "csr_matrix")
            self.assertEqual(profile["matrix"]["nnz"], 3)
            self.assertAlmostEqual(profile["matrix"]["density"], 0.5)

    def test_profiles_from_curation_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5ad_path = root / "dataset.h5ad"
            _write_dense_fixture(h5ad_path)
            manifest = {
                "schema_version": "1.0",
                "datasets": [
                    {
                        "dataset": "fixture",
                        "output": str(h5ad_path),
                        "output_sha256": "abc",
                    }
                ],
            }
            manifest_path = root / "curation-manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            output_dir = root / "profiles"
            index = profile_from_curation_manifest(manifest_path, output_dir)
            self.assertEqual(index["profile_count"], 1)
            self.assertTrue((output_dir / "fixture.profile.json").exists())
            self.assertTrue((output_dir / "profile-index.json").exists())

    def test_profiles_from_release_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5ad_path = root / "dataset.h5ad"
            _write_dense_fixture(h5ad_path)
            release = {
                "schema_version": "1.0",
                "release_id": "fixture_v0.1",
                "datasets": [
                    {
                        "dataset_id": "fixture",
                        "analysis_input": {
                            "path": str(h5ad_path),
                            "storage": "h5ad",
                            "mode": "passthrough",
                            "sha256": "abc",
                        },
                    }
                ],
            }
            release_path = root / "release.json"
            release_path.write_text(json.dumps(release), encoding="utf-8")
            output_dir = root / "profiles"
            index = profile_from_release_manifest(release_path, output_dir)
            self.assertEqual(index["release_id"], "fixture_v0.1")
            self.assertEqual(index["profile_count"], 1)


if __name__ == "__main__":
    unittest.main()
