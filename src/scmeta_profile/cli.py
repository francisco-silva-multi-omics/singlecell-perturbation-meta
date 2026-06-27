from __future__ import annotations

import argparse
from pathlib import Path
import sys

from scmeta_profile.profiler import (
    ProfileError,
    profile_from_curation_manifest,
    profile_from_release_manifest,
    profile_h5ad,
    write_json,
)


def _profile(args: argparse.Namespace) -> int:
    if args.release_manifest:
        index = profile_from_release_manifest(
            args.release_manifest,
            args.output_dir,
            compute_sha256=args.compute_sha256,
        )
        print(f"Wrote {index['profile_count']} profiles to {args.output_dir}")
        print(f"Wrote profile index to {args.output_dir / 'profile-index.json'}")
        return 0

    if args.curation_manifest:
        index = profile_from_curation_manifest(
            args.curation_manifest,
            args.output_dir,
            compute_sha256=args.compute_sha256,
        )
        print(f"Wrote {index['profile_count']} profiles to {args.output_dir}")
        print(f"Wrote profile index to {args.output_dir / 'profile-index.json'}")
        return 0

    if not args.input or not args.dataset_id or not args.output:
        raise ProfileError(
            "Use --curation-manifest/--output-dir or --input/--dataset-id/--output"
        )
    profile = profile_h5ad(
        args.input,
        dataset_id=args.dataset_id,
        compute_sha256=args.compute_sha256,
    )
    write_json(profile, args.output)
    print(f"Wrote profile to {args.output}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-profile")
    subparsers = parser.add_subparsers(dest="command", required=True)

    profile = subparsers.add_parser("profile", help="Profile H5AD files")
    profile.add_argument("--release-manifest", type=Path)
    profile.add_argument("--curation-manifest", type=Path)
    profile.add_argument("--output-dir", type=Path, default=Path("metadata/profiles"))
    profile.add_argument("--input", type=Path)
    profile.add_argument("--dataset-id")
    profile.add_argument("--output", type=Path)
    profile.add_argument(
        "--compute-sha256",
        action="store_true",
        help="Hash the H5AD file while profiling. Existing manifest hashes are reused by default.",
    )
    profile.set_defaults(handler=_profile)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (ProfileError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
