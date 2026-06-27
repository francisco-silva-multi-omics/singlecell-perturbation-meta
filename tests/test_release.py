from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from scmeta_release.build import build_release_manifest


def _write_h5ad(path: Path):
    with h5py.File(path, "w") as h5ad:
        h5ad.create_dataset("X", data=np.asarray([[1, 0], [0, 2]]))
        obs = h5ad.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        obs.create_dataset("cell_id", data=np.asarray(["c1", "c2"], dtype=h5py.string_dtype()))
        obs.create_dataset("perturbation", data=np.asarray(["control", "gene"], dtype=h5py.string_dtype()))
        var = h5ad.create_group("var")
        var.attrs["_index"] = "gene_id"
        var.create_dataset("gene_id", data=np.asarray(["g1", "g2"], dtype=h5py.string_dtype()))


class ReleaseTests(unittest.TestCase):
    def test_builds_release_manifest_from_curation_and_profiles(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5ad_path = root / "dataset.h5ad"
            _write_h5ad(h5ad_path)
            curation_path = root / "curation.json"
            curation_path.write_text(
                json.dumps(
                    {
                        "schema_version": "1.0",
                        "datasets": [
                            {
                                "dataset": "fixture",
                                "input": str(h5ad_path),
                                "output": str(h5ad_path),
                                "mode": "passthrough",
                                "cells_before": None,
                                "cells_after": None,
                                "cells_excluded": None,
                                "output_size_bytes": h5ad_path.stat().st_size,
                                "output_sha256": "sha",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            profile_path = root / "fixture.profile.json"
            profile_path.write_text("{}", encoding="utf-8")
            profile_index_path = root / "profile-index.json"
            profile_index_path.write_text(
                json.dumps(
                    {
                        "schema_version": "1.0",
                        "profiles": [
                            {
                                "dataset_id": "fixture",
                                "profile_path": str(profile_path),
                                "h5ad_path": str(h5ad_path),
                                "n_obs": 2,
                                "n_vars": 2,
                                "warnings": ["warning"],
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            output_path = root / "release.json"
            release = build_release_manifest(
                release_id="test_v0.1",
                curation_manifest_path=curation_path,
                profile_index_path=profile_index_path,
                output_path=output_path,
            )
            self.assertEqual(release["dataset_count"], 1)
            self.assertEqual(release["datasets"][0]["analysis_input"]["path"], str(h5ad_path))
            self.assertEqual(release["datasets"][0]["cell_counts"]["profile_n_obs"], 2)
            self.assertTrue(output_path.exists())


if __name__ == "__main__":
    unittest.main()
