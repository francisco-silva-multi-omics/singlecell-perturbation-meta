from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class ReleaseError(RuntimeError):
    """Raised when release-manifest inputs are inconsistent."""


ADMISSION_BY_DATASET = {
    "adamson": "pass_after_curation",
    "dixit": "pass_after_curation",
    "norman": "pass_documented_subset",
    "replogle": "pass_primary_matrix_verified",
    "sciplex3": "pass_after_curation",
}


def _read_json(path: Path) -> dict[str, Any] | list[Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _profile_by_dataset(profile_index_path: Path) -> dict[str, dict[str, Any]]:
    payload = _read_json(profile_index_path)
    if not isinstance(payload, dict) or "profiles" not in payload:
        raise ReleaseError("Profile index must contain a profiles list")
    result = {}
    for profile in payload["profiles"]:
        dataset_id = str(profile["dataset_id"])
        result[dataset_id] = profile
    return result


def _download_by_name(receipt_path: Path | None) -> dict[str, dict[str, Any]]:
    if receipt_path is None:
        return {}
    payload = _read_json(receipt_path)
    if not isinstance(payload, list):
        raise ReleaseError("Download receipt must be a JSON list")
    return {str(item["name"]): item for item in payload}


def _load_optional_json(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    payload = _read_json(path)
    if not isinstance(payload, dict):
        raise ReleaseError(f"Expected JSON object: {path}")
    return payload


def build_release_manifest(
    *,
    release_id: str,
    curation_manifest_path: Path,
    profile_index_path: Path,
    output_path: Path,
    download_receipt_path: Path | None = None,
    primary_audit_path: Path | None = None,
    replogle_verification_path: Path | None = None,
) -> dict[str, Any]:
    curation = _read_json(curation_manifest_path)
    if not isinstance(curation, dict) or "datasets" not in curation:
        raise ReleaseError("Curation manifest must contain a datasets list")

    profiles = _profile_by_dataset(profile_index_path)
    downloads = _download_by_name(download_receipt_path)
    primary_audit = primary_audit_path if primary_audit_path else None
    replogle_verification = _load_optional_json(replogle_verification_path)

    datasets = []
    for record in curation["datasets"]:
        dataset_id = str(record["dataset"])
        if dataset_id not in profiles:
            raise ReleaseError(f"Missing profile for dataset: {dataset_id}")
        final_path = _resolve_path(record["output"])
        if not final_path.exists():
            raise ReleaseError(f"Final analysis input does not exist: {final_path}")
        profile = profiles[dataset_id]
        profile_path = _resolve_path(profile["profile_path"])
        if not profile_path.exists():
            raise ReleaseError(f"Profile JSON does not exist: {profile_path}")

        source_name = Path(record["input"]).name
        download = downloads.get(source_name, {})
        file_size = final_path.stat().st_size
        expected_size = record.get("output_size_bytes")
        if expected_size is not None and int(expected_size) != file_size:
            raise ReleaseError(
                f"Size mismatch for {dataset_id}: manifest {expected_size}, file {file_size}"
            )

        dataset = {
            "dataset_id": dataset_id,
            "release_id": release_id,
            "analysis_input": {
                "path": str(final_path),
                "storage": "h5ad",
                "mode": record.get("mode"),
                "size_bytes": file_size,
                "sha256": record.get("output_sha256"),
            },
            "lineage": {
                "landing_input": str(_resolve_path(record["input"])),
                "curation_manifest": str(curation_manifest_path),
                "curation_mode": record.get("mode"),
                "profile": str(profile_path),
                "download_receipt": str(download_receipt_path) if download_receipt_path else None,
                "primary_audit": str(primary_audit) if primary_audit else None,
            },
            "cell_counts": {
                "before_curation": record.get("cells_before"),
                "after_curation": record.get("cells_after"),
                "excluded_by_curation": record.get("cells_excluded"),
                "profile_n_obs": profile.get("n_obs"),
                "profile_n_vars": profile.get("n_vars"),
            },
            "admission_status": ADMISSION_BY_DATASET.get(dataset_id, "unknown"),
            "warnings": profile.get("warnings", []),
        }
        if download:
            dataset["lineage"]["download_artifact"] = {
                "name": download.get("name"),
                "size_bytes": download.get("size_bytes"),
                "sha256": download.get("sha256"),
                "status": download.get("status"),
            }
        if dataset_id == "replogle" and replogle_verification:
            dataset["lineage"]["primary_matrix_verification"] = {
                "path": str(replogle_verification_path),
                "status": replogle_verification.get("status"),
                "matrix_values_match": replogle_verification.get("matrix_values_match"),
                "primary_matrix_sha256": replogle_verification.get("primary_matrix_sha256"),
                "harmonized_matrix_sha256": replogle_verification.get(
                    "harmonized_matrix_sha256"
                ),
            }
        datasets.append(dataset)

    release = {
        "schema_version": "1.0",
        "release_id": release_id,
        "release_type": "analysis_input_manifest",
        "source": {
            "curation_manifest": str(curation_manifest_path),
            "profile_index": str(profile_index_path),
            "download_receipt": str(download_receipt_path) if download_receipt_path else None,
            "primary_audit": str(primary_audit) if primary_audit else None,
            "replogle_verification": str(replogle_verification_path)
            if replogle_verification_path
            else None,
        },
        "dataset_count": len(datasets),
        "datasets": sorted(datasets, key=lambda item: item["dataset_id"]),
    }
    _write_json(release, output_path)
    return release
