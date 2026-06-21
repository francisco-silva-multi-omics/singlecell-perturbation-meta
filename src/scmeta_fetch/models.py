from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "1.0"


class ManifestError(ValueError):
    """Raised when a source manifest violates the local data contract."""


@dataclass(frozen=True)
class Checksum:
    algorithm: str
    value: str

    def __post_init__(self) -> None:
        if self.algorithm not in {"md5", "sha256"} or not self.value:
            raise ManifestError(
                f"Unsupported checksum: {self.algorithm}:{self.value}"
            )

    @classmethod
    def parse(cls, raw: str) -> "Checksum":
        try:
            algorithm, value = raw.lower().split(":", 1)
        except ValueError as exc:
            raise ManifestError(f"Invalid checksum: {raw!r}") from exc
        if algorithm not in {"md5", "sha256"} or not value:
            raise ManifestError(f"Unsupported checksum: {raw!r}")
        return cls(algorithm=algorithm, value=value)


@dataclass(frozen=True)
class Artifact:
    name: str
    url: str
    size_bytes: int
    upstream_checksum: Checksum

    def __post_init__(self) -> None:
        if not self.name or Path(self.name).name != self.name:
            raise ManifestError(f"Artifact name must be a basename: {self.name!r}")
        if not self.url.startswith("https://"):
            raise ManifestError(f"Artifact URL must use HTTPS: {self.url!r}")
        if self.size_bytes < 0:
            raise ManifestError(f"Artifact size cannot be negative: {self.name}")


@dataclass(frozen=True)
class DatasetManifest:
    schema_version: str
    generated_at: str
    source: dict[str, Any]
    release: dict[str, Any]
    artifacts: tuple[Artifact, ...]

    @classmethod
    def create(
        cls,
        *,
        source: dict[str, Any],
        release: dict[str, Any],
        artifacts: list[Artifact],
    ) -> "DatasetManifest":
        return cls(
            schema_version=SCHEMA_VERSION,
            generated_at=datetime.now(timezone.utc).isoformat(),
            source=source,
            release=release,
            artifacts=tuple(sorted(artifacts, key=lambda item: item.name)),
        )

    @classmethod
    def from_dict(cls, raw: dict[str, Any]) -> "DatasetManifest":
        if raw.get("schema_version") != SCHEMA_VERSION:
            raise ManifestError(
                f"Unsupported manifest schema: {raw.get('schema_version')!r}"
            )
        artifacts = []
        for item in raw.get("artifacts", []):
            checksum = item.get("upstream_checksum", {})
            artifacts.append(
                Artifact(
                    name=item["name"],
                    url=item["url"],
                    size_bytes=int(item["size_bytes"]),
                    upstream_checksum=Checksum(
                        algorithm=checksum["algorithm"], value=checksum["value"]
                    ),
                )
            )
        if not artifacts:
            raise ManifestError("Manifest contains no artifacts")
        names = [item.name for item in artifacts]
        if len(names) != len(set(names)):
            raise ManifestError("Manifest contains duplicate artifact names")
        return cls(
            schema_version=raw["schema_version"],
            generated_at=raw["generated_at"],
            source=dict(raw["source"]),
            release=dict(raw["release"]),
            artifacts=tuple(artifacts),
        )

    @classmethod
    def read(cls, path: Path) -> "DatasetManifest":
        with path.open("r", encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))

    def write(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = asdict(self)
        temporary = path.with_suffix(path.suffix + ".partial")
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        temporary.replace(path)
