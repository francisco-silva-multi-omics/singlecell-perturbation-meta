from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_release.build import ReleaseError, build_release_manifest


def _build(args: argparse.Namespace) -> int:
    release = build_release_manifest(
        release_id=args.release_id,
        curation_manifest_path=args.curation_manifest,
        profile_index_path=args.profile_index,
        download_receipt_path=args.download_receipt,
        primary_audit_path=args.primary_audit,
        replogle_verification_path=args.replogle_verification,
        output_path=args.output,
    )
    print(f"Wrote release {release['release_id']} with {release['dataset_count']} datasets to {args.output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-release")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="Build an analysis-input release manifest")
    build.add_argument("--release-id", required=True)
    build.add_argument("--curation-manifest", type=Path, required=True)
    build.add_argument("--profile-index", type=Path, required=True)
    build.add_argument("--download-receipt", type=Path)
    build.add_argument("--primary-audit", type=Path)
    build.add_argument("--replogle-verification", type=Path)
    build.add_argument("--output", type=Path, required=True)
    build.set_defaults(handler=_build)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (ReleaseError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
