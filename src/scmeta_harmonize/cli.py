from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_harmonize.harmonize import HarmonizeError, harmonize_release


def _harmonize(args: argparse.Namespace) -> int:
    manifest = harmonize_release(
        release_manifest_path=args.release_manifest,
        output_dir=args.output_dir,
    )
    print(
        f"Wrote harmonized metadata for {len(manifest['datasets'])} datasets to {args.output_dir}"
    )
    print(
        f"Wrote {manifest['interventions']['rows']} interventions and "
        f"{manifest['conditions']['rows']} conditions"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-harmonize")
    subparsers = parser.add_subparsers(dest="command", required=True)

    harmonize = subparsers.add_parser(
        "harmonize", help="Export canonical cell, gene, intervention, and condition metadata"
    )
    harmonize.add_argument("--release-manifest", type=Path, required=True)
    harmonize.add_argument("--output-dir", type=Path, required=True)
    harmonize.set_defaults(handler=_harmonize)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (HarmonizeError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
