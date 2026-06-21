from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import sys

from scmeta_fetch.download import DownloadError, download_artifact, write_receipt
from scmeta_fetch.models import DatasetManifest, ManifestError
from scmeta_fetch.zenodo import resolve_record


def _size_label(size_bytes: int) -> str:
    value = float(size_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024
    raise RuntimeError("unreachable")


def _read_selection(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, list) or not all(isinstance(item, str) for item in payload):
        raise ManifestError("Selection file must be a JSON array of artifact names")
    return set(payload)


def _catalog(args: argparse.Namespace) -> int:
    manifest = resolve_record(args.record_id)
    manifest.write(args.output)
    total = sum(item.size_bytes for item in manifest.artifacts)
    print(f"Wrote {len(manifest.artifacts)} artifacts ({_size_label(total)}) to {args.output}")
    return 0


def _selected_artifacts(args: argparse.Namespace, manifest: DatasetManifest):
    requested = set(args.include or [])
    if args.selection_file:
        requested.update(_read_selection(args.selection_file))
    if args.all:
        requested.update(item.name for item in manifest.artifacts)
    if not requested:
        raise ManifestError("Select files with --include/--selection-file, or explicitly pass --all")

    by_name = {item.name: item for item in manifest.artifacts}
    missing = requested - by_name.keys()
    if missing:
        raise ManifestError(f"Artifacts not present in manifest: {', '.join(sorted(missing))}")
    return [by_name[name] for name in sorted(requested)]


def _fetch(args: argparse.Namespace) -> int:
    manifest = DatasetManifest.read(args.manifest)
    artifacts = _selected_artifacts(args, manifest)
    total = sum(item.size_bytes for item in artifacts)
    print(f"Selected {len(artifacts)} artifacts totaling {_size_label(total)}")
    for artifact in artifacts:
        print(f"  {artifact.name}: {_size_label(artifact.size_bytes)}")
    if args.dry_run:
        return 0
    if not args.yes:
        raise ManifestError("Refusing to download without --yes; use --dry-run to inspect first")

    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                download_artifact,
                artifact,
                args.output_dir,
                timeout=args.timeout,
                retries=args.retries,
            ): artifact
            for artifact in artifacts
        }
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{result.status}: {result.name} ({result.sha256[:12]}...)")

    receipt = args.receipt or args.output_dir / "download-receipt.json"
    write_receipt(results, receipt)
    print(f"Wrote receipt to {receipt}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="scmeta-fetch")
    subparsers = parser.add_subparsers(dest="command", required=True)

    catalog = subparsers.add_parser("catalog", help="Resolve a source release to a manifest")
    catalog_subparsers = catalog.add_subparsers(dest="source", required=True)
    zenodo = catalog_subparsers.add_parser("zenodo", help="Resolve a Zenodo record")
    zenodo.add_argument("--record-id", required=True)
    zenodo.add_argument("--output", type=Path, required=True)
    zenodo.set_defaults(handler=_catalog)

    fetch = subparsers.add_parser("fetch", help="Download selected manifest artifacts")
    fetch.add_argument("--manifest", type=Path, required=True)
    fetch.add_argument("--output-dir", type=Path, required=True)
    fetch.add_argument("--include", action="append")
    fetch.add_argument("--selection-file", type=Path)
    fetch.add_argument("--all", action="store_true")
    fetch.add_argument("--workers", type=int, default=2, choices=range(1, 9))
    fetch.add_argument("--timeout", type=float, default=60)
    fetch.add_argument("--retries", type=int, default=3)
    fetch.add_argument("--receipt", type=Path)
    fetch.add_argument("--dry-run", action="store_true")
    fetch.add_argument("--yes", action="store_true")
    fetch.set_defaults(handler=_fetch)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (DownloadError, ManifestError, OSError, ValueError) as exc:
        parser.error(str(exc))
    return 2


if __name__ == "__main__":
    sys.exit(main())
