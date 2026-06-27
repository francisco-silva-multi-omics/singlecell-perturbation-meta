from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_aggregate.aggregate import aggregate_harmonized
from scmeta_harmonize.harmonize import harmonize_release


def _string_dataset(parent, name, values):
    dataset = parent.create_dataset(
        name, data=np.asarray(values, dtype=h5py.string_dtype("utf-8"))
    )
    dataset.attrs["encoding-type"] = "string-array"
    dataset.attrs["encoding-version"] = "0.2.0"


def _write_fixture(path: Path):
    with h5py.File(path, "w") as h5ad:
        h5ad.create_dataset(
            "X",
            data=np.asarray(
                [
                    [1, 0, 2],
                    [0, 3, 0],
                    [4, 0, 1],
                    [0, 2, 2],
                ],
                dtype=np.float32,
            ),
        )
        obs = h5ad.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        _string_dataset(obs, "cell_id", ["c1", "c2", "c3", "c4"])
        _string_dataset(obs, "cell_line", ["K562", "K562", "K562", "K562"])
        _string_dataset(obs, "perturbation", ["control", "geneA", "geneA", "control"])
        _string_dataset(obs, "perturbation_type", ["CRISPR", "CRISPR", "CRISPR", "CRISPR"])
        _string_dataset(obs, "target", ["", "geneA", "geneA", ""])
        _string_dataset(obs, "guide_id", ["ctrl", "g1", "g1", "ctrl"])
        var = h5ad.create_group("var")
        var.attrs["_index"] = "gene_symbol"
        _string_dataset(var, "gene_symbol", ["g1", "g2", "g3"])


class AggregateTests(unittest.TestCase):
    def test_aggregates_dense_counts_by_intervention(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5ad_path = root / "fixture.h5ad"
            _write_fixture(h5ad_path)
            release_path = root / "release.json"
            release_path.write_text(
                json.dumps(
                    {
                        "schema_version": "1.0",
                        "release_id": "fixture_v0.1",
                        "datasets": [
                            {
                                "dataset_id": "fixture",
                                "analysis_input": {"path": str(h5ad_path)},
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            harmonized_dir = root / "harmonized"
            harmonize_manifest = harmonize_release(
                release_manifest_path=release_path,
                output_dir=harmonized_dir,
            )
            aggregate_dir = root / "aggregate"
            aggregate_manifest = aggregate_harmonized(
                harmonization_manifest_path=harmonized_dir / "harmonization-manifest.json",
                output_dir=aggregate_dir,
                group_by="intervention",
                chunk_size=2,
            )
            self.assertEqual(len(aggregate_manifest["datasets"]), 1)
            matrix_path = aggregate_dir / "intervention" / "fixture" / "pseudobulk.h5"
            with h5py.File(matrix_path, "r") as handle:
                data = handle["X/data"][()]
                indices = handle["X/indices"][()]
                indptr = handle["X/indptr"][()]
                dense = np.zeros(tuple(handle["X"].attrs["shape"]))
                for row in range(dense.shape[0]):
                    start, end = indptr[row], indptr[row + 1]
                    dense[row, indices[start:end]] = data[start:end]
                group_ids = [value.decode() for value in handle["group_id"][()]]
            with (aggregate_dir / "intervention" / "fixture" / "groups.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                groups = list(csv.DictReader(handle))
            by_perturbation = {
                row["perturbation_id"]: dense[group_ids.index(row["group_id"])]
                for row in groups
            }
            np.testing.assert_array_equal(by_perturbation["control"], [1, 2, 4])
            np.testing.assert_array_equal(by_perturbation["geneA"], [4, 3, 1])


if __name__ == "__main__":
    unittest.main()
