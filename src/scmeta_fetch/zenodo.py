from __future__ import annotations

from typing import Any, Callable

from scmeta_fetch.http import get_json
from scmeta_fetch.models import Artifact, Checksum, DatasetManifest, ManifestError


ZENODO_API = "https://zenodo.org/api/records"


def resolve_record(
    record_id: str,
    *,
    fetch_json: Callable[[str], dict[str, Any]] = get_json,
) -> DatasetManifest:
    if not record_id.isdigit():
        raise ManifestError(f"Zenodo record ID must be numeric: {record_id!r}")

    api_url = f"{ZENODO_API}/{record_id}"
    record = fetch_json(api_url)
    files = record.get("files") or []
    artifacts = [
        Artifact(
            name=item["key"],
            url=item["links"]["self"],
            size_bytes=int(item["size"]),
            upstream_checksum=Checksum.parse(item["checksum"]),
        )
        for item in files
    ]
    if not artifacts:
        raise ManifestError(f"Zenodo record {record_id} contains no downloadable files")

    metadata = record.get("metadata", {})
    links = record.get("links", {})
    license_data = metadata.get("license") or {}
    return DatasetManifest.create(
        source={
            "name": "zenodo",
            "record_id": str(record.get("id", record_id)),
            "api_url": links.get("self", api_url),
            "landing_page": links.get("self_html", f"https://zenodo.org/records/{record_id}"),
        },
        release={
            "title": metadata.get("title"),
            "doi": record.get("doi"),
            "version": metadata.get("version"),
            "license": license_data.get("id"),
            "modified": record.get("modified"),
        },
        artifacts=artifacts,
    )

