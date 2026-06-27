from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_aggregate.aggregate import aggregate_harmonized
from scmeta_effects.effects import compute_effects, compute_pooled_effects
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
                    [10, 0],
                    [10, 0],
                    [0, 10],
                    [0, 10],
                ],
                dtype=np.float32,
            ),
        )
        obs = h5ad.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        _string_dataset(obs, "cell_id", ["c1", "c2", "c3", "c4"])
        _string_dataset(obs, "perturbation", ["control", "control", "geneA", "geneA"])
        _string_dataset(obs, "perturbation_type", ["CRISPR", "CRISPR", "CRISPR", "CRISPR"])
        _string_dataset(obs, "target", ["", "", "geneA", "geneA"])
        _string_dataset(obs, "guide_id", ["ctrl", "ctrl", "g1", "g1"])
        var = h5ad.create_group("var")
        var.attrs["_index"] = "gene_symbol"
        _string_dataset(var, "gene_symbol", ["g1", "g2"])


class EffectTests(unittest.TestCase):
    def test_computes_treated_minus_control_effects(self):
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
            harmonize_release(release_manifest_path=release_path, output_dir=harmonized_dir)
            aggregate_dir = root / "aggregate"
            aggregate_harmonized(
                harmonization_manifest_path=harmonized_dir / "harmonization-manifest.json",
                output_dir=aggregate_dir,
                group_by="intervention",
                chunk_size=2,
            )
            effect_dir = root / "effects"
            manifest = compute_effects(
                pseudobulk_manifest_path=aggregate_dir / "intervention" / "pseudobulk-manifest.json",
                output_dir=effect_dir,
                min_cells=1,
                top_genes=2,
            )
            self.assertEqual(manifest["datasets"][0]["contrasts"], 1)
            with h5py.File(effect_dir / "fixture" / "effects.h5", "r") as handle:
                delta = handle["delta_log1p_cpm"][()]
                features = [value.decode() for value in handle["feature_id"][()]]
            self.assertEqual(features, ["g1", "g2"])
            self.assertLess(delta[0, 0], 0)
            self.assertGreater(delta[0, 1], 0)
            with (effect_dir / "fixture" / "contrasts.csv").open(newline="", encoding="utf-8") as handle:
                contrasts = list(csv.DictReader(handle))
            self.assertEqual(contrasts[0]["perturbation_id"], "geneA")
            self.assertTrue((effect_dir / "fixture" / "top_genes.csv.gz").exists())

    def test_computes_pooled_condition_effects(self):
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
            harmonize_release(release_manifest_path=release_path, output_dir=harmonized_dir)
            aggregate_dir = root / "aggregate"
            aggregate_harmonized(
                harmonization_manifest_path=harmonized_dir / "harmonization-manifest.json",
                output_dir=aggregate_dir,
                group_by="condition",
                chunk_size=2,
            )

            with (harmonized_dir / "conditions.csv").open(newline="", encoding="utf-8") as handle:
                conditions = list(csv.DictReader(handle))
            control_id = next(row["condition_id"] for row in conditions if row["is_control"] == "true")
            treated_id = next(row["condition_id"] for row in conditions if row["is_control"] != "true")
            pooled_path = root / "pooled-contrasts.csv"
            with pooled_path.open("w", newline="", encoding="utf-8") as handle:
                fieldnames = [
                    "contrast_id",
                    "dataset_id",
                    "perturbation_id",
                    "dose",
                    "dose_unit",
                    "time",
                    "time_unit",
                    "cell_line",
                    "cell_type",
                    "disease",
                    "status",
                    "reason",
                    "match_rule",
                    "treated_condition_ids",
                    "control_condition_ids",
                    "treated_cell_count",
                    "control_cell_count",
                    "treated_condition_count",
                    "control_condition_count",
                    "matched_strata_count",
                ]
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(
                    {
                        "contrast_id": "pooled_contrast_00000001",
                        "dataset_id": "fixture",
                        "perturbation_id": "geneA",
                        "status": "pass",
                        "treated_condition_ids": treated_id,
                        "control_condition_ids": control_id,
                        "treated_cell_count": "2",
                        "control_cell_count": "2",
                        "treated_condition_count": "1",
                        "control_condition_count": "1",
                        "matched_strata_count": "1",
                    }
                )

            effect_dir = root / "pooled-effects"
            manifest = compute_pooled_effects(
                pseudobulk_manifest_path=aggregate_dir / "condition" / "pseudobulk-manifest.json",
                pooled_contrasts_path=pooled_path,
                output_dir=effect_dir,
                top_genes=2,
            )
            self.assertEqual(manifest["datasets"][0]["contrasts"], 1)
            with h5py.File(effect_dir / "fixture" / "effects.h5", "r") as handle:
                delta = handle["delta_log1p_cpm"][()]
                features = [value.decode() for value in handle["feature_id"][()]]
            self.assertEqual(features, ["g1", "g2"])
            self.assertLess(delta[0, 0], 0)
            self.assertGreater(delta[0, 1], 0)
            with (effect_dir / "fixture" / "contrasts.csv").open(newline="", encoding="utf-8") as handle:
                contrasts = list(csv.DictReader(handle))
            self.assertEqual(contrasts[0]["treated_condition_ids"], treated_id)
            self.assertEqual(contrasts[0]["control_condition_ids"], control_id)


if __name__ == "__main__":
    unittest.main()
