from __future__ import annotations

import hashlib
from io import BytesIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scmeta_fetch.download import DownloadError, download_artifact
from scmeta_fetch.models import Artifact, Checksum


class FakeResponse:
    def __init__(self, content: bytes, status: int) -> None:
        self._stream = BytesIO(content)
        self.status = status

    def read(self, size: int = -1) -> bytes:
        return self._stream.read(size)

    def getcode(self) -> int:
        return self.status

    def __enter__(self) -> "FakeResponse":
        return self

    def __exit__(self, *_args) -> None:
        return None


def make_artifact(content: bytes) -> Artifact:
    return Artifact(
        name="example.bin",
        url="https://example.test/example.bin",
        size_bytes=len(content),
        upstream_checksum=Checksum("md5", hashlib.md5(content).hexdigest()),
    )


class DownloadTests(unittest.TestCase):
    def test_downloads_verifies_and_promotes_atomically(self) -> None:
        content = b"single-cell-data"
        artifact = make_artifact(content)
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir)
            with patch(
                "scmeta_fetch.download.open_with_retries",
                return_value=FakeResponse(content, 200),
            ):
                result = download_artifact(artifact, output)

            self.assertEqual((output / artifact.name).read_bytes(), content)
            self.assertFalse((output / f"{artifact.name}.partial").exists())
            self.assertEqual(result.status, "downloaded")
            self.assertEqual(result.sha256, hashlib.sha256(content).hexdigest())

    def test_resumes_a_partial_download(self) -> None:
        content = b"abcdefghij"
        artifact = make_artifact(content)
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir)
            (output / f"{artifact.name}.partial").write_bytes(content[:4])
            with patch(
                "scmeta_fetch.download.open_with_retries",
                return_value=FakeResponse(content[4:], 206),
            ) as opener:
                result = download_artifact(artifact, output)

            request = opener.call_args.args[0]
            self.assertEqual(request.headers["Range"], "bytes=4-")
            self.assertEqual((output / artifact.name).read_bytes(), content)
            self.assertEqual(result.status, "resumed")

    def test_keeps_bad_content_quarantined_as_partial(self) -> None:
        expected = b"good-content"
        artifact = make_artifact(expected)
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir)
            with patch(
                "scmeta_fetch.download.open_with_retries",
                return_value=FakeResponse(b"bad-content!", 200),
            ):
                with self.assertRaises(DownloadError):
                    download_artifact(artifact, output)

            self.assertFalse((output / artifact.name).exists())
            self.assertTrue((output / f"{artifact.name}.partial").exists())


if __name__ == "__main__":
    unittest.main()

