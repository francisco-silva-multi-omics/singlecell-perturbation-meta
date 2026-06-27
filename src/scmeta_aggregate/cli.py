from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_aggregate.aggregate import AggregateError, aggregate_harmonized


def _aggregate(args: argparse.Namespace) -> int:
    include = set(args.include_dataset or []) or None
    manifest = aggregate_harmonized(
        harmonization_manifest_path=args.harmonization_manifest,
        output_dir=args.output_dir,
        group_by=args.group_by,
        include_datasets=include,
        chunk_size=args.chunk_size,
    )
    print(
        f"Wrote {manifest['group_by']} pseudobulk for "
        f"{len(manifest['datasets'])} datasets to {args.output_dir / manifest['group_by']}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-aggregate")
    subparsers = parser.add_subparsers(dest="command", required=True)

    aggregate = subparsers.add_parser("aggregate", help="Build pseudobulk count matrices")
    aggregate.add_argument("--harmonization-manifest", type=Path, required=True)
    aggregate.add_argument("--output-dir", type=Path, required=True)
    aggregate.add_argument(
        "--group-by",
        choices=("intervention", "bio-condition", "condition"),
        default="intervention",
    )
    aggregate.add_argument("--include-dataset", action="append")
    aggregate.add_argument("--chunk-size", type=int, default=4096)
    aggregate.set_defaults(handler=_aggregate)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (AggregateError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
