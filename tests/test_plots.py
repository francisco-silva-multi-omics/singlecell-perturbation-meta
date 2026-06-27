from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_plots.plots import plot_mvp_summary


class PlotTests(unittest.TestCase):
    def test_writes_svg_summary_figures(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            release = {
                "datasets": [
                    {
                        "dataset_id": "fixture",
                        "cell_counts": {"profile_n_obs": 4, "profile_n_vars": 3},
                    }
                ]
            }
            release_path = root / "release.json"
            release_path.write_text(json.dumps(release), encoding="utf-8")
            qc_path = root / "qc.json"
            qc_path.write_text(
                json.dumps(
                    {
                        "dataset_status_counts": {"warn": 1},
                        "perturbation_status_counts": {"pass": 1, "control": 1},
                    }
                ),
                encoding="utf-8",
            )
            effects_dir = root / "effects" / "fixture"
            effects_dir.mkdir(parents=True)
            effects_h5 = effects_dir / "effects.h5"
            with h5py.File(effects_h5, "w") as handle:
                handle.create_dataset("delta_log1p_cpm", data=np.asarray([[1.0, -2.0, 0.5]]))
            top_path = effects_dir / "top_genes.csv.gz"
            with gzip.open(top_path, "wt", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["gene_id", "effect_size"])
                writer.writeheader()
                writer.writerow({"gene_id": "geneA", "effect_size": "2.0"})
            contrasts_path = effects_dir / "contrasts.csv"
            with contrasts_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["contrast_id", "status"])
                writer.writeheader()
                writer.writerow({"contrast_id": "c1", "status": "pass"})
            effects_path = root / "effects.json"
            effects_path.write_text(
                json.dumps(
                    {
                        "model": "treated_minus_matched_control_delta",
                        "datasets": [
                            {
                                "dataset_id": "fixture",
                                "contrasts": 1,
                                "controls": 1,
                                "matrix": {"path": str(effects_h5)},
                                "contrast_table": {"path": str(contrasts_path)},
                                "top_genes": {"path": str(top_path)},
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            output_dir = root / "figures"
            manifest = plot_mvp_summary(
                release_manifest_path=release_path,
                qc_summary_path=qc_path,
                effects_manifest_path=effects_path,
                output_dir=output_dir,
            )
            self.assertEqual(len(manifest["figures"]), 3)
            self.assertEqual(manifest["effects_model"], "treated_minus_matched_control_delta")
            self.assertEqual(manifest["datasets"][0]["status_counts"]["pass"], 1)
            for path in manifest["figures"].values():
                self.assertTrue(Path(path).exists())
                self.assertIn("<svg", Path(path).read_text(encoding="utf-8"))
            self.assertTrue((output_dir / "visual_summary.md").exists())


if __name__ == "__main__":
    unittest.main()
