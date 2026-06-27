from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request

from scmeta_fetch.http import USER_AGENT, open_with_retries
from scmeta_fetch.models import Artifact


class DownloadError(RuntimeError):
    """Raised when a download or integrity check fails."""


@dataclass(frozen=True)
class DownloadResult:
    name: str
    path: str
    size_bytes: int
    sha256: str
    status: str


def _digests(path: Path, upstream_algorithm: str) -> tuple[str, str]:
    sha256 = hashlib.sha256()
    upstream = hashlib.new(upstream_algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            sha256.update(chunk)
            upstream.update(chunk)
    return sha256.hexdigest(), upstream.hexdigest()


def _verify(path: Path, artifact: Artifact) -> str:
    actual_size = path.stat().st_size
    if actual_size != artifact.size_bytes:
        raise DownloadError(
            f"Size mismatch for {artifact.name}: expected {artifact.size_bytes}, got {actual_size}"
        )
    sha256, upstream = _digests(path, artifact.upstream_checksum.algorithm)
    if upstream != artifact.upstream_checksum.value:
        raise DownloadError(
            f"Checksum mismatch for {artifact.name}: expected "
            f"{artifact.upstream_checksum.algorithm}:{artifact.upstream_checksum.value}, "
            f"got {upstream}"
        )
    return sha256


def download_artifact(
    artifact: Artifact,
    output_dir: Path,
    *,
    timeout: float = 60,
    retries: int = 3,
) -> DownloadResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / artifact.name
    partial = output_dir / f"{artifact.name}.partial"

    if destination.exists():
        sha256 = _verify(destination, artifact)
        return DownloadResult(
            artifact.name, str(destination), artifact.size_bytes, sha256, "already_verified"
        )

    offset = partial.stat().st_size if partial.exists() else 0
    if offset > artifact.size_bytes:
        raise DownloadError(f"Partial file is larger than expected: {partial}")

    # Zenodo's content endpoint returns 406 when an explicit octet-stream
    # Accept header is sent, even though the successful response uses that type.
    headers = {"User-Agent": USER_AGENT}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = Request(artifact.url, headers=headers)

    try:
        response = open_with_retries(request, timeout=timeout, retries=retries)
    except HTTPError as exc:
        if exc.code == 416 and offset == artifact.size_bytes:
            sha256 = _verify(partial, artifact)
            os.replace(partial, destination)
            return DownloadResult(
                artifact.name, str(destination), artifact.size_bytes, sha256, "resumed"
            )
        raise

    with response:
        response_status = getattr(response, "status", response.getcode())
        append = offset > 0 and response_status == 206
        mode = "ab" if append else "wb"
        with partial.open(mode) as handle:
            while chunk := response.read(8 * 1024 * 1024):
                handle.write(chunk)

    sha256 = _verify(partial, artifact)
    os.replace(partial, destination)
    return DownloadResult(
        artifact.name,
        str(destination),
        artifact.size_bytes,
        sha256,
        "resumed" if offset and append else "downloaded",
    )


def write_receipt(results: list[DownloadResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump([asdict(item) for item in sorted(results, key=lambda x: x.name)], handle, indent=2)
        handle.write("\n")
    temporary.replace(path)
