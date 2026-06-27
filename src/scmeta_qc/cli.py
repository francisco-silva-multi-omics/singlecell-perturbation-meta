from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_qc.qc import QCError, run_qc


def _qc(args: argparse.Namespace) -> int:
    summary = run_qc(
        release_manifest_path=args.release_manifest,
        output_dir=args.output_dir,
        min_cells=args.min_cells,
        min_control_cells=args.min_control_cells,
    )
    print(
        f"Wrote QC for {summary['dataset_count']} datasets and "
        f"{summary['perturbation_count']} perturbations to {args.output_dir}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-qc")
    subparsers = parser.add_subparsers(dest="command", required=True)

    qc = subparsers.add_parser("qc", help="Run dataset and perturbation QC gates")
    qc.add_argument("--release-manifest", type=Path, required=True)
    qc.add_argument("--output-dir", type=Path, required=True)
    qc.add_argument("--min-cells", type=int, default=50)
    qc.add_argument("--min-control-cells", type=int, default=50)
    qc.set_defaults(handler=_qc)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (QCError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
