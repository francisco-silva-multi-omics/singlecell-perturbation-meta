from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_effects.effects import EffectError, compute_effects, compute_pooled_effects


def _effects(args: argparse.Namespace) -> int:
    include = set(args.include_dataset or []) or None
    manifest = compute_effects(
        pseudobulk_manifest_path=args.pseudobulk_manifest,
        output_dir=args.output_dir,
        include_datasets=include,
        scale=args.scale,
        min_cells=args.min_cells,
        top_genes=args.top_genes,
    )
    total = sum(item["contrasts"] for item in manifest["datasets"])
    print(f"Wrote effects for {len(manifest['datasets'])} datasets and {total} contrasts to {args.output_dir}")
    return 0


def _pooled(args: argparse.Namespace) -> int:
    include = set(args.include_dataset or []) or None
    statuses = set(args.include_status or []) or None
    manifest = compute_pooled_effects(
        pseudobulk_manifest_path=args.pseudobulk_manifest,
        pooled_contrasts_path=args.pooled_contrasts,
        output_dir=args.output_dir,
        include_statuses=statuses,
        include_datasets=include,
        scale=args.scale,
        top_genes=args.top_genes,
    )
    total = sum(item["contrasts"] for item in manifest["datasets"])
    print(
        f"Wrote pooled effects for {len(manifest['datasets'])} datasets "
        f"and {total} contrasts to {args.output_dir}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-effects")
    subparsers = parser.add_subparsers(dest="command", required=True)

    effects = subparsers.add_parser("effects", help="Compute treated-vs-control effect matrices")
    effects.add_argument("--pseudobulk-manifest", type=Path, required=True)
    effects.add_argument("--output-dir", type=Path, required=True)
    effects.add_argument("--include-dataset", action="append")
    effects.add_argument("--scale", type=float, default=1_000_000.0)
    effects.add_argument("--min-cells", type=int, default=50)
    effects.add_argument("--top-genes", type=int, default=100)
    effects.set_defaults(handler=_effects)

    pooled = subparsers.add_parser(
        "pooled", help="Compute effects from pooled condition-level contrast definitions"
    )
    pooled.add_argument("--pseudobulk-manifest", type=Path, required=True)
    pooled.add_argument("--pooled-contrasts", type=Path, required=True)
    pooled.add_argument("--output-dir", type=Path, required=True)
    pooled.add_argument("--include-dataset", action="append")
    pooled.add_argument(
        "--include-status",
        action="append",
        help="Contrast status to include; repeatable. Defaults to pass and warn.",
    )
    pooled.add_argument("--scale", type=float, default=1_000_000.0)
    pooled.add_argument("--top-genes", type=int, default=100)
    pooled.set_defaults(handler=_pooled)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (EffectError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
