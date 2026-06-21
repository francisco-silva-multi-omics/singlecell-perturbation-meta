from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from scmeta_fetch.models import Artifact, Checksum, DatasetManifest, ManifestError


class DatasetManifestTests(unittest.TestCase):
    def test_round_trip(self) -> None:
        manifest = DatasetManifest.create(
            source={"name": "test"},
            release={"version": "1"},
            artifacts=[
                Artifact("b.h5ad", "https://example.test/b", 2, Checksum("md5", "abc")),
                Artifact("a.h5ad", "https://example.test/a", 1, Checksum("md5", "def")),
            ],
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "manifest.json"
            manifest.write(path)
            loaded = DatasetManifest.read(path)

        self.assertEqual([item.name for item in loaded.artifacts], ["a.h5ad", "b.h5ad"])
        self.assertEqual(loaded.release["version"], "1")

    def test_rejects_duplicate_names(self) -> None:
        payload = {
            "schema_version": "1.0",
            "generated_at": "2026-01-01T00:00:00+00:00",
            "source": {},
            "release": {},
            "artifacts": [
                {
                    "name": "same.h5ad",
                    "url": "https://example.test/one",
                    "size_bytes": 1,
                    "upstream_checksum": {"algorithm": "md5", "value": "a"},
                },
                {
                    "name": "same.h5ad",
                    "url": "https://example.test/two",
                    "size_bytes": 1,
                    "upstream_checksum": {"algorithm": "md5", "value": "b"},
                },
            ],
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "manifest.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(ManifestError):
                DatasetManifest.read(path)

    def test_rejects_unsupported_checksum(self) -> None:
        with self.assertRaises(ManifestError):
            Checksum("crc32", "1234")


if __name__ == "__main__":
    unittest.main()
