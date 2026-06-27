from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_meta.meta import MetaError, summarize_pooled_effects


def _summarize(args: argparse.Namespace) -> int:
    statuses = set(args.include_status or []) or None
    manifest = summarize_pooled_effects(
        effects_manifest_path=args.effects_manifest,
        output_dir=args.output_dir,
        include_statuses=statuses,
        top_n=args.top_n,
    )
    print(
        f"Wrote meta summaries for {manifest['gene_count']} genes "
        f"and {manifest['included_top_gene_hits']} top-effect hits to {args.output_dir}"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-meta")
    subparsers = parser.add_subparsers(dest="command", required=True)

    summarize = subparsers.add_parser(
        "summarize", help="Summarize pooled effect recurrence and sign consistency"
    )
    summarize.add_argument("--effects-manifest", type=Path, required=True)
    summarize.add_argument("--output-dir", type=Path, required=True)
    summarize.add_argument(
        "--include-status",
        action="append",
        help="Contrast status to include; repeatable. Defaults to pass and warn.",
    )
    summarize.add_argument("--top-n", type=int, default=100)
    summarize.set_defaults(handler=_summarize)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (MetaError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
