from __future__ import annotations

import unittest

from scmeta_fetch.zenodo import resolve_record


class ZenodoTests(unittest.TestCase):
    def test_resolves_record_to_local_contract(self) -> None:
        response = {
            "id": 123,
            "doi": "10.5281/zenodo.123",
            "modified": "2026-01-01T00:00:00Z",
            "metadata": {
                "title": "Example",
                "version": "1.2",
                "license": {"id": "cc-by-4.0"},
            },
            "links": {
                "self": "https://zenodo.org/api/records/123",
                "self_html": "https://zenodo.org/records/123",
            },
            "files": [
                {
                    "key": "example.h5ad",
                    "size": 42,
                    "checksum": "md5:abc123",
                    "links": {"self": "https://zenodo.org/api/files/example/content"},
                }
            ],
        }

        manifest = resolve_record("123", fetch_json=lambda _: response)

        self.assertEqual(manifest.source["record_id"], "123")
        self.assertEqual(manifest.release["version"], "1.2")
        self.assertEqual(manifest.artifacts[0].name, "example.h5ad")
        self.assertEqual(manifest.artifacts[0].upstream_checksum.value, "abc123")


if __name__ == "__main__":
    unittest.main()

