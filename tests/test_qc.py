from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_profile.profiler import profile_h5ad, write_json
from scmeta_qc.qc import run_qc


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
                    [0, 1, 0],
                    [2, 0, 0],
                    [0, 3, 1],
                    [1, 1, 0],
                ]
            ),
        )
        obs = h5ad.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        _string_dataset(obs, "cell_id", [f"c{i}" for i in range(5)])
        _string_dataset(obs, "perturbation", ["control", "control", "geneA", "geneA", "geneB"])
        _string_dataset(obs, "control_status", ["negative_control", "negative_control", "targeting", "targeting", "targeting"])
        _string_dataset(obs, "batch", ["b1", "b2", "b1", "b2", "b1"])
        var = h5ad.create_group("var")
        var.attrs["_index"] = "gene_id"
        _string_dataset(var, "gene_id", ["g1", "g2", "g3"])


class QCTests(unittest.TestCase):
    def test_runs_dataset_and_perturbation_qc(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5ad_path = root / "fixture.h5ad"
            _write_fixture(h5ad_path)
            profile = profile_h5ad(h5ad_path, dataset_id="fixture")
            profile_path = root / "fixture.profile.json"
            write_json(profile, profile_path)
            release = {
                "schema_version": "1.0",
                "release_id": "fixture_v0.1",
                "datasets": [
                    {
                        "dataset_id": "fixture",
                        "admission_status": "pass_after_curation",
                        "analysis_input": {
                            "path": str(h5ad_path),
                            "storage": "h5ad",
                            "mode": "curated",
                            "sha256": "sha",
                        },
                        "lineage": {"profile": str(profile_path)},
                    }
                ],
            }
            release_path = root / "release.json"
            release_path.write_text(json.dumps(release), encoding="utf-8")
            output_dir = root / "qc"
            summary = run_qc(
                release_manifest_path=release_path,
                output_dir=output_dir,
                min_cells=2,
                min_control_cells=2,
            )
            self.assertEqual(summary["dataset_count"], 1)
            self.assertEqual(summary["perturbation_count"], 3)
            dataset_qc = json.loads((output_dir / "dataset-qc.json").read_text())
            self.assertEqual(dataset_qc["datasets"][0]["control_cells"], 2)
            perturbation_qc = json.loads((output_dir / "perturbation-qc.json").read_text())
            by_perturbation = {
                row["perturbation_id"]: row for row in perturbation_qc["perturbations"]
            }
            self.assertEqual(by_perturbation["control"]["status"], "control")
            self.assertEqual(by_perturbation["geneA"]["status"], "pass")
            self.assertEqual(by_perturbation["geneB"]["status"], "fail")
            self.assertTrue((output_dir / "dataset-qc.csv").exists())
            self.assertTrue((output_dir / "perturbation-qc.csv").exists())


if __name__ == "__main__":
    unittest.main()
