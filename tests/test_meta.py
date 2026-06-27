from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from scmeta_meta.meta import summarize_pooled_effects


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_gzip_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


class MetaTests(unittest.TestCase):
    def test_summarizes_recurrence_and_sign_consistency(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            datasets = []
            top_fields = [
                "dataset_id",
                "contrast_id",
                "perturbation_id",
                "gene_id",
                "effect_size",
                "treated_log1p_cpm",
                "control_log1p_cpm",
                "rank_abs_effect",
            ]
            contrast_fields = ["contrast_id", "status"]
            fixtures = {
                "ds1": [
                    {"contrast_id": "c1", "status": "pass", "gene_id": "geneA", "effect_size": "2.0"},
                    {"contrast_id": "c1", "status": "pass", "gene_id": "geneB", "effect_size": "-1.0"},
                    {"contrast_id": "c2", "status": "fail", "gene_id": "geneA", "effect_size": "5.0"},
                ],
                "ds2": [
                    {"contrast_id": "c3", "status": "warn", "gene_id": "geneA", "effect_size": "1.5"},
                    {"contrast_id": "c3", "status": "warn", "gene_id": "geneC", "effect_size": "-3.0"},
                ],
            }
            for dataset_id, rows in fixtures.items():
                dataset_dir = root / dataset_id
                contrasts = []
                top_rows = []
                for row in rows:
                    contrasts.append(
                        {"contrast_id": row["contrast_id"], "status": row["status"]}
                    )
                    top_rows.append(
                        {
                            "dataset_id": dataset_id,
                            "contrast_id": row["contrast_id"],
                            "perturbation_id": "pert",
                            "gene_id": row["gene_id"],
                            "effect_size": row["effect_size"],
                            "treated_log1p_cpm": "",
                            "control_log1p_cpm": "",
                            "rank_abs_effect": "1",
                        }
                    )
                contrast_path = dataset_dir / "contrasts.csv"
                top_path = dataset_dir / "top_genes.csv.gz"
                _write_csv(contrast_path, contrasts, contrast_fields)
                _write_gzip_csv(top_path, top_rows, top_fields)
                datasets.append(
                    {
                        "dataset_id": dataset_id,
                        "contrasts": len({row["contrast_id"] for row in rows}),
                        "contrast_table": {"path": str(contrast_path)},
                        "top_genes": {"path": str(top_path)},
                    }
                )
            manifest_path = root / "effects-manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "release_id": "fixture",
                        "model": "treated_minus_matched_control_delta",
                        "datasets": datasets,
                    }
                ),
                encoding="utf-8",
            )
            output_dir = root / "meta"
            manifest = summarize_pooled_effects(
                effects_manifest_path=manifest_path,
                output_dir=output_dir,
                top_n=10,
            )
            self.assertEqual(manifest["gene_count"], 3)
            self.assertEqual(manifest["included_top_gene_hits"], 4)
            with (output_dir / "gene-consistency.csv").open(newline="", encoding="utf-8") as handle:
                rows = {row["gene_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(rows["geneA"]["dataset_count"], "2")
            self.assertEqual(rows["geneA"]["pass_hits"], "1")
            self.assertEqual(rows["geneA"]["warn_hits"], "1")
            self.assertEqual(float(rows["geneA"]["sign_consistency"]), 1.0)
            self.assertNotIn("5.0", rows["geneA"]["mean_effect"])
            self.assertTrue((output_dir / "meta-summary.csv").exists())
            self.assertTrue((output_dir / "meta-manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
