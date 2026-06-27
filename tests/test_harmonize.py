from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_harmonize.harmonize import harmonize_release


def _string_dataset(parent, name, values):
    dataset = parent.create_dataset(
        name, data=np.asarray(values, dtype=h5py.string_dtype("utf-8"))
    )
    dataset.attrs["encoding-type"] = "string-array"
    dataset.attrs["encoding-version"] = "0.2.0"


def _write_fixture(path: Path):
    with h5py.File(path, "w") as h5ad:
        h5ad.create_dataset("X", data=np.asarray([[1, 0], [0, 1], [2, 2]]))
        obs = h5ad.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        _string_dataset(obs, "cell_id", ["c1", "c2", "c3"])
        _string_dataset(obs, "cell_line", ["K562", "K562", "K562"])
        _string_dataset(obs, "celltype", ["lymphoblasts", "lymphoblasts", "lymphoblasts"])
        _string_dataset(obs, "disease", ["leukemia", "leukemia", "leukemia"])
        _string_dataset(obs, "organism", ["human", "human", "human"])
        _string_dataset(obs, "perturbation", ["control", "geneA", "geneA"])
        _string_dataset(obs, "perturbation_type", ["CRISPR", "CRISPR", "CRISPR"])
        _string_dataset(obs, "target", ["", "geneA", "geneA"])
        _string_dataset(obs, "guide_id", ["ctrl", "g1", "g1"])
        _string_dataset(obs, "batch", ["b1", "b1", "b2"])
        obs.create_dataset("ncounts", data=np.asarray([10, 20, 30]))
        obs.create_dataset("ngenes", data=np.asarray([2, 2, 2]))
        obs.create_dataset("percent_mito", data=np.asarray([1.0, 2.0, 3.0]))
        obs.create_dataset("percent_ribo", data=np.asarray([4.0, 5.0, 6.0]))

        var = h5ad.create_group("var")
        var.attrs["_index"] = "gene_symbol"
        _string_dataset(var, "gene_symbol", ["A", "B"])
        _string_dataset(var, "ensembl_id", ["ENSG1", "ENSG2"])
        _string_dataset(var, "chr", ["chr1", "chr2"])
        var.create_dataset("start", data=np.asarray([1, 2]))
        var.create_dataset("end", data=np.asarray([10, 20]))


def _read_gzip_csv(path: Path):
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class HarmonizeTests(unittest.TestCase):
    def test_harmonizes_release_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5ad_path = root / "fixture.h5ad"
            _write_fixture(h5ad_path)
            release = {
                "schema_version": "1.0",
                "release_id": "fixture_v0.1",
                "datasets": [
                    {
                        "dataset_id": "fixture",
                        "analysis_input": {"path": str(h5ad_path)},
                    }
                ],
            }
            release_path = root / "release.json"
            release_path.write_text(json.dumps(release), encoding="utf-8")
            output_dir = root / "harmonized"
            manifest = harmonize_release(
                release_manifest_path=release_path,
                output_dir=output_dir,
            )

            self.assertEqual(manifest["release_id"], "fixture_v0.1")
            obs_rows = _read_gzip_csv(output_dir / "obs" / "fixture.csv.gz")
            self.assertEqual(len(obs_rows), 3)
            self.assertEqual(obs_rows[0]["is_control"], "true")
            self.assertEqual(obs_rows[1]["perturbation_target"], "geneA")
            var_rows = _read_gzip_csv(output_dir / "var" / "fixture.csv.gz")
            self.assertEqual(var_rows[0]["ensembl_id"], "ENSG1")
            with (output_dir / "conditions.csv").open(newline="", encoding="utf-8") as handle:
                conditions = list(csv.DictReader(handle))
            with (output_dir / "interventions.csv").open(newline="", encoding="utf-8") as handle:
                interventions = list(csv.DictReader(handle))
            self.assertEqual(len(interventions), 2)
            self.assertEqual(sum(int(row["cell_count"]) for row in conditions), 3)
            self.assertTrue((output_dir / "harmonization-manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
