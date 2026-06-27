from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np


class ProfileError(RuntimeError):
    """Raised when an H5AD file cannot be profiled with the expected contract."""


OBS_REQUIRED = {
    "perturbation",
}

OBS_RECOMMENDED = {
    "batch",
    "batch_id",
    "cell_type",
    "control_status",
    "dose",
    "drug",
    "drug_dose",
    "drug_time",
    "guide_id",
    "n_counts",
    "n_genes",
    "nperts",
    "perturbation_target",
    "perturbation_type",
    "replicate",
    "time",
}

VAR_RECOMMENDED = {
    "ensembl_id",
    "gene_id",
    "gene_name",
    "gene_symbol",
}

PERTURBATION_COLUMNS = (
    "perturbation",
    "perturbation_target",
    "guide_id",
    "drug",
    "drug_name",
    "condition",
)

CONTROL_COLUMNS = ("control_status", "is_control", "perturbation")
BATCH_COLUMNS = ("batch", "batch_id", "gem_group", "replicate", "replicate_id", "plate")
DOSE_TIME_COLUMNS = ("dose", "drug_dose", "time", "drug_time", "timepoint")


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)!r}")


def write_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _decode_scalar(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.bytes_):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def _decode_array(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in "SO":
        return np.asarray([_decode_scalar(value) for value in values], dtype=object)
    return values


def _attr(group: h5py.Group | h5py.Dataset, key: str) -> Any:
    if key not in group.attrs:
        return None
    return _decode_scalar(group.attrs[key])


def _attrs(group: h5py.Group | h5py.Dataset) -> dict[str, Any]:
    return {key: _decode_scalar(value) for key, value in group.attrs.items()}


def _dataframe_index_name(group: h5py.Group) -> str | None:
    value = _attr(group, "_index")
    return None if value is None else str(value)


def _column_order(group: h5py.Group) -> list[str]:
    values = group.attrs.get("column-order", [])
    return [str(_decode_scalar(value)) for value in values]


def _read_column(group: h5py.Group, name: str) -> np.ndarray:
    obj = group[name]
    if isinstance(obj, h5py.Dataset):
        return _decode_array(obj[()])
    if _attr(obj, "encoding-type") != "categorical":
        raise ProfileError(f"Unsupported dataframe column encoding for {name}")
    categories = _read_column(obj, "categories")
    codes = obj["codes"][()]
    return np.asarray(
        [None if code < 0 else _decode_scalar(categories[int(code)]) for code in codes],
        dtype=object,
    )


def _count_missing(values: np.ndarray) -> int:
    if values.dtype.kind in "fc":
        return int(np.isnan(values).sum())
    missing = 0
    for value in values:
        if value is None:
            missing += 1
            continue
        text = str(value)
        if text == "" or text.lower() in {"nan", "none", "null", "na"}:
            missing += 1
    return missing


def _value_counts(values: np.ndarray, *, limit: int = 20) -> list[dict[str, Any]]:
    counter: Counter[str] = Counter()
    for value in values:
        if value is None:
            key = "<missing>"
        else:
            key = str(_decode_scalar(value))
            if key == "" or key.lower() in {"nan", "none", "null", "na"}:
                key = "<missing>"
        counter[key] += 1
    return [
        {"value": value, "count": count}
        for value, count in counter.most_common(limit)
    ]


def _column_profile(group: h5py.Group, name: str) -> dict[str, Any]:
    obj = group[name]
    base: dict[str, Any] = {"name": name}
    if isinstance(obj, h5py.Dataset):
        base.update(
            {
                "storage": "dataset",
                "dtype": str(obj.dtype),
                "shape": list(obj.shape),
                "encoding_type": _attr(obj, "encoding-type"),
            }
        )
        if obj.ndim == 1:
            values = _decode_array(obj[()])
            base["missing_count"] = _count_missing(values)
            base["unique_count"] = int(len(set(map(str, values.tolist()))))
        return base

    encoding = _attr(obj, "encoding-type")
    base.update({"storage": "group", "encoding_type": encoding})
    if encoding == "categorical" and "codes" in obj and "categories" in obj:
        codes = obj["codes"][()]
        categories = _read_column(obj, "categories")
        base.update(
            {
                "dtype": str(obj["codes"].dtype),
                "shape": list(codes.shape),
                "category_count": int(len(categories)),
                "missing_count": int((codes < 0).sum()),
                "categories_sample": [
                    str(_decode_scalar(value)) for value in categories[:20]
                ],
            }
        )
    return base


def _dataframe_profile(group: h5py.Group, *, summary_columns: Iterable[str]) -> dict[str, Any]:
    columns = sorted(group.keys())
    profile: dict[str, Any] = {
        "encoding_type": _attr(group, "encoding-type"),
        "encoding_version": _attr(group, "encoding-version"),
        "index": _dataframe_index_name(group),
        "column_order": _column_order(group),
        "columns": [_column_profile(group, name) for name in columns],
        "column_names": columns,
        "summaries": {},
    }

    if profile["index"] and profile["index"] in group:
        profile["n_rows"] = int(len(_read_column(group, profile["index"])))
    else:
        profile["n_rows"] = None

    for name in summary_columns:
        if name not in group:
            continue
        values = _read_column(group, name)
        profile["summaries"][name] = {
            "missing_count": _count_missing(values),
            "unique_count": int(len(set(map(str, values.tolist())))),
            "top_values": _value_counts(values),
        }
    return profile


def _matrix_shape(obj: h5py.Dataset | h5py.Group) -> tuple[int, int]:
    if isinstance(obj, h5py.Dataset):
        if obj.ndim != 2:
            raise ProfileError(f"X dataset must be 2-dimensional, found {obj.ndim}")
        return int(obj.shape[0]), int(obj.shape[1])
    if "shape" not in obj.attrs:
        raise ProfileError("Sparse X group is missing shape attribute")
    shape = tuple(int(value) for value in obj.attrs["shape"])
    if len(shape) != 2:
        raise ProfileError(f"X shape must have two dimensions, found {shape}")
    return shape


def _is_integer_like(values: np.ndarray) -> bool:
    if values.size == 0:
        return True
    if values.dtype.kind in "iu":
        return True
    if values.dtype.kind not in "f":
        return False
    return bool(np.all(np.isfinite(values)) and np.allclose(values, np.round(values)))


def _dense_sample(dataset: h5py.Dataset, *, max_rows: int = 64) -> np.ndarray:
    row_count = dataset.shape[0]
    if row_count == 0:
        return np.asarray([], dtype=dataset.dtype)
    take = min(max_rows, row_count)
    if take == row_count:
        rows = np.arange(row_count)
    else:
        rows = np.unique(np.linspace(0, row_count - 1, num=take, dtype=np.int64))
    blocks = [np.asarray(dataset[int(row) : int(row) + 1]) for row in rows]
    return np.concatenate(blocks, axis=0)


def _sparse_data_sample(data: h5py.Dataset, *, max_values: int = 200_000) -> np.ndarray:
    if len(data) <= max_values:
        return np.asarray(data[()])
    window = max(1, max_values // 4)
    starts = [
        0,
        max(0, len(data) // 3 - window // 2),
        max(0, (2 * len(data)) // 3 - window // 2),
        max(0, len(data) - window),
    ]
    samples = [np.asarray(data[start : min(start + window, len(data))]) for start in starts]
    return np.concatenate(samples)


def _matrix_profile(h5ad: h5py.File) -> dict[str, Any]:
    if "X" not in h5ad:
        raise ProfileError("H5AD is missing X")
    x = h5ad["X"]
    shape = _matrix_shape(x)
    profile: dict[str, Any] = {
        "shape": list(shape),
        "n_obs": shape[0],
        "n_vars": shape[1],
    }
    if isinstance(x, h5py.Dataset):
        sample = _dense_sample(x)
        nonzero_sample = int(np.count_nonzero(sample))
        profile.update(
            {
                "encoding": "dense",
                "dtype": str(x.dtype),
                "sampled_values": int(sample.size),
                "sampled_nonzero": nonzero_sample,
                "sampled_density": None
                if sample.size == 0
                else nonzero_sample / float(sample.size),
                "integer_like_sample": _is_integer_like(sample),
                "nnz": None,
                "density": None,
            }
        )
        return profile

    encoding = _attr(x, "encoding-type")
    if "data" not in x:
        raise ProfileError("Sparse X group is missing data")
    data = x["data"]
    nnz = int(len(data))
    sample = _sparse_data_sample(data)
    total = shape[0] * shape[1]
    profile.update(
        {
            "encoding": encoding,
            "dtype": str(data.dtype),
            "nnz": nnz,
            "density": None if total == 0 else nnz / float(total),
            "sampled_values": int(sample.size),
            "integer_like_sample": _is_integer_like(sample),
        }
    )
    return profile


def _counts_locations(h5ad: h5py.File) -> dict[str, Any]:
    layers = []
    if "layers" in h5ad and isinstance(h5ad["layers"], h5py.Group):
        layers = sorted(h5ad["layers"].keys())
    return {
        "x_present": "X" in h5ad,
        "raw_present": "raw" in h5ad,
        "layers": layers,
        "counts_location_inferred": "X",
    }


def _warnings(profile: dict[str, Any]) -> list[str]:
    warnings = []
    obs_columns = set(profile["obs"]["column_names"])
    var_columns = set(profile["var"]["column_names"])
    if profile["obs"].get("index"):
        obs_columns.add(profile["obs"]["index"])
    if profile["var"].get("index"):
        var_columns.add(profile["var"]["index"])
    missing_required = sorted(OBS_REQUIRED - obs_columns)
    if missing_required:
        warnings.append(f"missing required obs columns: {', '.join(missing_required)}")
    missing_recommended = sorted(OBS_RECOMMENDED - obs_columns)
    if missing_recommended:
        warnings.append(
            "missing recommended obs columns: " + ", ".join(missing_recommended)
        )
    if not (VAR_RECOMMENDED & var_columns):
        warnings.append(
            "missing recommended gene identifier columns: "
            + ", ".join(sorted(VAR_RECOMMENDED))
        )
    if not profile["counts"]["raw_present"] and not profile["counts"]["layers"]:
        warnings.append("no raw object or count layer found; X is treated as counts")
    if not profile["matrix"]["integer_like_sample"]:
        warnings.append("sampled X values are not integer-like")
    obs_rows = profile["obs"].get("n_rows")
    var_rows = profile["var"].get("n_rows")
    if obs_rows is not None and obs_rows != profile["matrix"]["n_obs"]:
        warnings.append("obs row count does not match X rows")
    if var_rows is not None and var_rows != profile["matrix"]["n_vars"]:
        warnings.append("var row count does not match X columns")
    return warnings


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def profile_h5ad(
    path: Path,
    *,
    dataset_id: str,
    source_record: dict[str, Any] | None = None,
    compute_sha256: bool = False,
) -> dict[str, Any]:
    if not path.exists():
        raise ProfileError(f"H5AD file does not exist: {path}")
    with h5py.File(path, "r") as h5ad:
        if "obs" not in h5ad or "var" not in h5ad:
            raise ProfileError("H5AD must contain obs and var groups")
        profile = {
            "schema_version": "1.0",
            "dataset_id": dataset_id,
            "path": str(path),
            "file_size_bytes": path.stat().st_size,
            "source_record": source_record or {},
            "h5ad_attrs": _attrs(h5ad),
            "matrix": _matrix_profile(h5ad),
            "counts": _counts_locations(h5ad),
            "obs": _dataframe_profile(
                h5ad["obs"],
                summary_columns=(
                    set(PERTURBATION_COLUMNS)
                    | set(CONTROL_COLUMNS)
                    | set(BATCH_COLUMNS)
                    | set(DOSE_TIME_COLUMNS)
                ),
            ),
            "var": _dataframe_profile(
                h5ad["var"],
                summary_columns=("gene_symbol", "gene_name", "ensembl_id"),
            ),
        }
    if compute_sha256:
        profile["file_sha256"] = _sha256(path)
    elif source_record:
        profile["file_sha256"] = source_record.get("output_sha256")
    profile["warnings"] = _warnings(profile)
    return profile


def _safe_stem(value: str) -> str:
    keep = []
    for char in value.lower():
        if char.isalnum():
            keep.append(char)
        elif char in {"-", "_", "."}:
            keep.append(char)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "dataset"


def _resolve_path(path: str | Path, *, base_dir: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    cwd_candidate = Path.cwd() / candidate
    if cwd_candidate.exists():
        return cwd_candidate
    return base_dir / candidate


def profile_from_curation_manifest(
    manifest_path: Path,
    output_dir: Path,
    *,
    compute_sha256: bool = False,
) -> dict[str, Any]:
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if "datasets" not in manifest or not isinstance(manifest["datasets"], list):
        raise ProfileError("Curation manifest must contain a datasets list")

    base_dir = manifest_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles = []
    for record in manifest["datasets"]:
        dataset_id = str(record["dataset"])
        source_path = _resolve_path(record["output"], base_dir=base_dir)
        profile = profile_h5ad(
            source_path,
            dataset_id=dataset_id,
            source_record=record,
            compute_sha256=compute_sha256,
        )
        profile_path = output_dir / f"{_safe_stem(dataset_id)}.profile.json"
        write_json(profile, profile_path)
        profiles.append(
            {
                "dataset_id": dataset_id,
                "profile_path": str(profile_path),
                "h5ad_path": str(source_path),
                "n_obs": profile["matrix"]["n_obs"],
                "n_vars": profile["matrix"]["n_vars"],
                "warnings": profile["warnings"],
            }
        )

    index = {
        "schema_version": "1.0",
        "source_manifest": str(manifest_path),
        "profile_count": len(profiles),
        "profiles": profiles,
    }
    write_json(index, output_dir / "profile-index.json")
    return index


def profile_from_release_manifest(
    release_manifest_path: Path,
    output_dir: Path,
    *,
    compute_sha256: bool = False,
) -> dict[str, Any]:
    with release_manifest_path.open("r", encoding="utf-8") as handle:
        release = json.load(handle)
    if "datasets" not in release or not isinstance(release["datasets"], list):
        raise ProfileError("Release manifest must contain a datasets list")

    output_dir.mkdir(parents=True, exist_ok=True)
    profiles = []
    for record in release["datasets"]:
        dataset_id = str(record["dataset_id"])
        source_path = _resolve_path(record["analysis_input"]["path"], base_dir=Path.cwd())
        profile = profile_h5ad(
            source_path,
            dataset_id=dataset_id,
            source_record=record,
            compute_sha256=compute_sha256,
        )
        profile_path = output_dir / f"{_safe_stem(dataset_id)}.profile.json"
        write_json(profile, profile_path)
        profiles.append(
            {
                "dataset_id": dataset_id,
                "profile_path": str(profile_path),
                "h5ad_path": str(source_path),
                "n_obs": profile["matrix"]["n_obs"],
                "n_vars": profile["matrix"]["n_vars"],
                "warnings": profile["warnings"],
            }
        )

    index = {
        "schema_version": "1.0",
        "source_manifest": str(release_manifest_path),
        "release_id": release.get("release_id"),
        "profile_count": len(profiles),
        "profiles": profiles,
    }
    write_json(index, output_dir / "profile-index.json")
    return index
