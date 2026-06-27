from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

from scmeta_curate.adamson import load_repair_metadata
from scmeta_curate.h5ad import CurationError, decode_array, filter_h5ad_rows
from scmeta_curate.replogle import verify_replogle, write_verification


FILES = {
    "adamson": "AdamsonWeissman2016_GSM2406681_10X010.h5ad",
    "dixit": "DixitRegev2016_K562_TFs_7_days.h5ad",
    "norman": "NormanWeissman2019_filtered.h5ad",
    "replogle": "ReplogleWeissman2022_K562_essential.h5ad",
    "sciplex3": "SrivatsanTrapnell2020_sciplex3.h5ad",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _filter_by_perturbation(input_path: Path, *, literal_nan: bool = False):
    with h5py.File(input_path, "r") as h5ad:
        perturbation = decode_array(h5ad["obs"], "perturbation")
    keep = np.asarray([value is not None for value in perturbation], dtype=bool)
    if literal_nan:
        keep &= np.asarray([str(value).lower() != "nan" for value in perturbation])
    nperts = np.asarray(
        [0 if str(value).lower() == "control" else 1 for value in perturbation],
        dtype=np.int64,
    )
    return keep, {"nperts": nperts}


def _curate(args: argparse.Namespace) -> int:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []

    adamson_input = args.input_dir / FILES["adamson"]
    keep, updates = load_repair_metadata(adamson_input, args.geo_cache)
    adamson_output = args.output_dir / FILES["adamson"]
    filter_h5ad_rows(adamson_input, adamson_output, keep, obs_updates=updates)
    records.append(("adamson", adamson_input, adamson_output, len(keep), int(keep.sum())))

    for key, literal_nan in (("dixit", True), ("sciplex3", False)):
        input_path = args.input_dir / FILES[key]
        output_path = args.output_dir / FILES[key]
        keep, updates = _filter_by_perturbation(input_path, literal_nan=literal_nan)
        filter_h5ad_rows(input_path, output_path, keep, obs_updates=updates)
        records.append((key, input_path, output_path, len(keep), int(keep.sum())))

    for key in ("norman", "replogle"):
        input_path = args.input_dir / FILES[key]
        records.append((key, input_path, input_path, None, None))

    manifest = {"schema_version": "1.0", "datasets": []}
    for key, input_path, output_path, before, after in records:
        manifest["datasets"].append(
            {
                "dataset": key,
                "input": str(input_path),
                "output": str(output_path),
                "mode": "curated" if input_path != output_path else "passthrough",
                "cells_before": before,
                "cells_after": after,
                "cells_excluded": None if before is None else before - after,
                "output_size_bytes": output_path.stat().st_size,
                "output_sha256": _sha256(output_path),
            }
        )
    manifest_path = args.output_dir / "curation-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote curation manifest: {manifest_path}")
    return 0


def _verify_replogle(args: argparse.Namespace) -> int:
    result = verify_replogle(args.primary, args.harmonized)
    write_verification(result, args.output)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "pass" else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-curate")
    subparsers = parser.add_subparsers(dest="command", required=True)
    curate = subparsers.add_parser("curate", help="Repair and filter the MVP H5AD cohort")
    curate.add_argument("--input-dir", type=Path, required=True)
    curate.add_argument("--output-dir", type=Path, required=True)
    curate.add_argument(
        "--geo-cache",
        type=Path,
        default=Path("data/raw/geo/GSE90546/GSM2406681"),
    )
    curate.set_defaults(handler=_curate)

    verify = subparsers.add_parser(
        "verify-replogle", help="Compare primary and harmonized Replogle H5AD files"
    )
    verify.add_argument("--primary", type=Path, required=True)
    verify.add_argument("--harmonized", type=Path, required=True)
    verify.add_argument("--output", type=Path, required=True)
    verify.set_defaults(handler=_verify_replogle)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (CurationError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2
