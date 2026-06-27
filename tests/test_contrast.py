from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from scmeta_contrast.contrast import build_contrasts, build_pooled_contrasts


class ContrastTests(unittest.TestCase):
    def test_builds_strict_and_fallback_contrasts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            conditions_path = root / "conditions.csv"
            fieldnames = [
                "condition_id",
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
                "cell_count",
            ]
            rows = [
                {
                    "condition_id": "c_control_b1",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "cell_type": "",
                    "disease": "leukemia",
                    "perturbation_id": "control",
                    "is_control": "true",
                    "batch_id": "b1",
                    "cell_count": "100",
                },
                {
                    "condition_id": "c_gene_b1",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "cell_type": "",
                    "disease": "leukemia",
                    "perturbation_id": "geneA",
                    "is_control": "false",
                    "batch_id": "b1",
                    "cell_count": "80",
                },
                {
                    "condition_id": "c_gene_b2",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "cell_type": "",
                    "disease": "leukemia",
                    "perturbation_id": "geneB",
                    "is_control": "false",
                    "batch_id": "b2",
                    "cell_count": "80",
                },
                {
                    "condition_id": "c_low",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "cell_type": "",
                    "disease": "leukemia",
                    "perturbation_id": "geneC",
                    "is_control": "false",
                    "batch_id": "b1",
                    "cell_count": "10",
                },
            ]
            with conditions_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                for row in rows:
                    writer.writerow({field: row.get(field, "") for field in fieldnames})
            manifest_path = root / "harmonization.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "release_id": "fixture",
                        "conditions": {"path": str(conditions_path)},
                    }
                ),
                encoding="utf-8",
            )
            output_dir = root / "contrasts"
            summary = build_contrasts(
                harmonization_manifest_path=manifest_path,
                output_dir=output_dir,
                min_treated_cells=50,
                min_control_cells=50,
            )
            self.assertEqual(summary["contrast_count"], 3)
            with (output_dir / "condition-contrasts.csv").open(newline="", encoding="utf-8") as handle:
                contrasts = {row["treated_condition_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(
                contrasts["c_gene_b1"]["match_rule"], "strict_context_batch_replicate"
            )
            self.assertEqual(contrasts["c_gene_b1"]["status"], "pass")
            self.assertEqual(contrasts["c_gene_b2"]["match_rule"], "bio_context")
            self.assertEqual(contrasts["c_gene_b2"]["status"], "pass")
            self.assertEqual(contrasts["c_low"]["status"], "fail")

    def test_builds_pooled_replicate_aware_contrasts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            conditions_path = root / "conditions.csv"
            fieldnames = [
                "condition_id",
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
                "cell_count",
            ]
            rows = [
                {
                    "condition_id": "ctrl_b1",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "disease": "leukemia",
                    "perturbation_id": "control",
                    "is_control": "true",
                    "batch_id": "b1",
                    "replicate_id": "r1",
                    "cell_count": "100",
                },
                {
                    "condition_id": "ctrl_b2",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "disease": "leukemia",
                    "perturbation_id": "control",
                    "is_control": "true",
                    "batch_id": "b2",
                    "replicate_id": "r2",
                    "cell_count": "100",
                },
                {
                    "condition_id": "gene_a_b1",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "disease": "leukemia",
                    "perturbation_id": "geneA",
                    "is_control": "false",
                    "batch_id": "b1",
                    "replicate_id": "r1",
                    "cell_count": "30",
                },
                {
                    "condition_id": "gene_a_b2",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "disease": "leukemia",
                    "perturbation_id": "geneA",
                    "is_control": "false",
                    "batch_id": "b2",
                    "replicate_id": "r2",
                    "cell_count": "30",
                },
                {
                    "condition_id": "gene_b_b3",
                    "dataset_id": "ds",
                    "cell_line": "K562",
                    "disease": "leukemia",
                    "perturbation_id": "geneB",
                    "is_control": "false",
                    "batch_id": "b3",
                    "replicate_id": "r3",
                    "cell_count": "70",
                },
            ]
            with conditions_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                for row in rows:
                    writer.writerow({field: row.get(field, "") for field in fieldnames})
            manifest_path = root / "harmonization.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "release_id": "fixture",
                        "conditions": {"path": str(conditions_path)},
                    }
                ),
                encoding="utf-8",
            )
            output_dir = root / "pooled"
            summary = build_pooled_contrasts(
                harmonization_manifest_path=manifest_path,
                output_dir=output_dir,
                min_treated_cells=50,
                min_control_cells=50,
                min_strata=2,
            )
            self.assertEqual(summary["contrast_count"], 2)
            with (output_dir / "pooled-contrasts.csv").open(newline="", encoding="utf-8") as handle:
                contrasts = {row["perturbation_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(contrasts["geneA"]["status"], "pass")
            self.assertEqual(contrasts["geneA"]["treated_cell_count"], "60")
            self.assertEqual(contrasts["geneA"]["matched_strata_count"], "2")
            self.assertEqual(contrasts["geneB"]["status"], "fail")
            self.assertIn("no matched control strata", contrasts["geneB"]["reason"])


if __name__ == "__main__":
    unittest.main()
