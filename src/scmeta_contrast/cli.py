from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_contrast.contrast import ContrastError, build_contrasts, build_pooled_contrasts


def _build(args: argparse.Namespace) -> int:
    summary = build_contrasts(
        harmonization_manifest_path=args.harmonization_manifest,
        output_dir=args.output_dir,
        min_treated_cells=args.min_treated_cells,
        min_control_cells=args.min_control_cells,
    )
    print(
        f"Wrote {summary['contrast_count']} condition contrasts to "
        f"{summary['output']}"
    )
    return 0


def _build_pooled(args: argparse.Namespace) -> int:
    summary = build_pooled_contrasts(
        harmonization_manifest_path=args.harmonization_manifest,
        output_dir=args.output_dir,
        min_treated_cells=args.min_treated_cells,
        min_control_cells=args.min_control_cells,
        min_strata=args.min_strata,
    )
    print(
        f"Wrote {summary['contrast_count']} pooled contrasts to "
        f"{summary['output']}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-contrast")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser(
        "build", help="Build explicit treated/control condition contrasts"
    )
    build.add_argument("--harmonization-manifest", type=Path, required=True)
    build.add_argument("--output-dir", type=Path, required=True)
    build.add_argument("--min-treated-cells", type=int, default=50)
    build.add_argument("--min-control-cells", type=int, default=50)
    build.set_defaults(handler=_build)

    pooled = subparsers.add_parser(
        "build-pooled", help="Build pooled replicate-aware treated/control contrasts"
    )
    pooled.add_argument("--harmonization-manifest", type=Path, required=True)
    pooled.add_argument("--output-dir", type=Path, required=True)
    pooled.add_argument("--min-treated-cells", type=int, default=50)
    pooled.add_argument("--min-control-cells", type=int, default=50)
    pooled.add_argument("--min-strata", type=int, default=2)
    pooled.set_defaults(handler=_build_pooled)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (ContrastError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
