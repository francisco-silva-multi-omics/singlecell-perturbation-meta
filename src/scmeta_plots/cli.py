from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_plots.plots import PlotError, plot_mvp_summary


def _plot(args: argparse.Namespace) -> int:
    manifest = plot_mvp_summary(
        release_manifest_path=args.release_manifest,
        qc_summary_path=args.qc_summary,
        effects_manifest_path=args.effects_manifest,
        output_dir=args.output_dir,
    )
    print(f"Wrote {len(manifest['figures'])} figures and summary report to {args.output_dir}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-plots")
    subparsers = parser.add_subparsers(dest="command", required=True)

    plot = subparsers.add_parser("mvp-summary", help="Plot MVP pipeline and effect summaries")
    plot.add_argument("--release-manifest", type=Path, required=True)
    plot.add_argument("--qc-summary", type=Path, required=True)
    plot.add_argument("--effects-manifest", type=Path, required=True)
    plot.add_argument("--output-dir", type=Path, required=True)
    plot.set_defaults(handler=_plot)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (PlotError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
