from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from typing import Any


class MetaError(RuntimeError):
    """Raised when meta-analysis inputs are missing or inconsistent."""


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise MetaError(f"Expected JSON object: {path}")
    return payload


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _write_csv(rows: list[dict[str, Any]], path: Path, fieldnames: list[str]) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    temporary.replace(path)
    return {"path": str(path), "rows": len(rows), "sha256": _sha256(path)}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _safe_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return 0.0


def _contrast_statuses(path: Path, include_statuses: set[str]) -> dict[str, str]:
    statuses = {}
    for row in _read_csv(path):
        status = row.get("status", "")
        if status in include_statuses:
            statuses[row["contrast_id"]] = status
    return statuses


def _iter_top_genes(path: Path):
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def _sign(value: float) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def _score(row: dict[str, Any]) -> float:
    return (
        float(row["dataset_count"])
        * float(row["sign_consistency"])
        * math.log10(float(row["hit_count"]) + 1.0)
        * float(row["mean_abs_effect"])
    )


def summarize_pooled_effects(
    *,
    effects_manifest_path: Path,
    output_dir: Path,
    include_statuses: set[str] | None = None,
    top_n: int = 100,
) -> dict[str, Any]:
    statuses = include_statuses or {"pass", "warn"}
    manifest = _read_json(effects_manifest_path)
    if "datasets" not in manifest or not isinstance(manifest["datasets"], list):
        raise MetaError("Effects manifest must contain a datasets list")

    gene_hits: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "effects": [],
            "datasets": set(),
            "contrasts": set(),
            "positive_hits": 0,
            "negative_hits": 0,
            "status_counts": Counter(),
            "dataset_counts": Counter(),
            "dataset_effects": defaultdict(list),
            "max_abs_effect": 0.0,
        }
    )
    dataset_rows = []
    total_input_hits = 0
    total_used_hits = 0

    for dataset in manifest["datasets"]:
        dataset_id = dataset["dataset_id"]
        contrast_path = _resolve_path(dataset["contrast_table"]["path"])
        top_path = _resolve_path(dataset["top_genes"]["path"])
        contrast_status = _contrast_statuses(contrast_path, statuses)
        used_hits = 0
        input_hits = 0
        for row in _iter_top_genes(top_path):
            input_hits += 1
            status = contrast_status.get(row["contrast_id"])
            if status is None:
                continue
            gene = row["gene_id"]
            effect = _safe_float(row.get("effect_size", "0"))
            hit = gene_hits[gene]
            hit["effects"].append(effect)
            hit["datasets"].add(dataset_id)
            hit["contrasts"].add(row["contrast_id"])
            hit["status_counts"][status] += 1
            hit["dataset_counts"][dataset_id] += 1
            hit["dataset_effects"][dataset_id].append(effect)
            hit["max_abs_effect"] = max(hit["max_abs_effect"], abs(effect))
            sign = _sign(effect)
            if sign > 0:
                hit["positive_hits"] += 1
            elif sign < 0:
                hit["negative_hits"] += 1
            used_hits += 1
        dataset_rows.append(
            {
                "dataset_id": dataset_id,
                "input_contrasts": dataset.get("contrasts", 0),
                "included_contrasts": len(contrast_status),
                "input_top_gene_hits": input_hits,
                "included_top_gene_hits": used_hits,
            }
        )
        total_input_hits += input_hits
        total_used_hits += used_hits

    gene_rows = []
    for gene, values in gene_hits.items():
        effects = values["effects"]
        hit_count = len(effects)
        if hit_count == 0:
            continue
        positive_hits = int(values["positive_hits"])
        negative_hits = int(values["negative_hits"])
        sign_consistency = max(positive_hits, negative_hits) / hit_count
        dataset_medians = {
            dataset_id: median(dataset_effects)
            for dataset_id, dataset_effects in values["dataset_effects"].items()
        }
        positive_datasets = sum(1 for value in dataset_medians.values() if value > 0)
        negative_datasets = sum(1 for value in dataset_medians.values() if value < 0)
        dataset_sign_consistency = (
            max(positive_datasets, negative_datasets) / len(dataset_medians)
            if dataset_medians
            else 0.0
        )
        row = {
            "gene_id": gene,
            "dataset_count": len(values["datasets"]),
            "contrast_count": len(values["contrasts"]),
            "hit_count": hit_count,
            "pass_hits": int(values["status_counts"].get("pass", 0)),
            "warn_hits": int(values["status_counts"].get("warn", 0)),
            "positive_hits": positive_hits,
            "negative_hits": negative_hits,
            "sign_consistency": sign_consistency,
            "dataset_sign_consistency": dataset_sign_consistency,
            "mean_effect": sum(effects) / hit_count,
            "median_effect": median(effects),
            "mean_abs_effect": sum(abs(value) for value in effects) / hit_count,
            "max_abs_effect": float(values["max_abs_effect"]),
            "datasets": ";".join(sorted(values["datasets"])),
        }
        row["conserved_response_score"] = _score(row)
        gene_rows.append(row)

    gene_rows.sort(
        key=lambda row: (
            row["dataset_count"],
            row["dataset_sign_consistency"],
            row["sign_consistency"],
            row["hit_count"],
            row["mean_abs_effect"],
        ),
        reverse=True,
    )
    summary_rows = []
    for rank, row in enumerate(gene_rows[:top_n], start=1):
        ranked = {"rank": rank}
        ranked.update(row)
        summary_rows.append(ranked)

    output_dir.mkdir(parents=True, exist_ok=True)
    gene_fieldnames = [
        "gene_id",
        "dataset_count",
        "contrast_count",
        "hit_count",
        "pass_hits",
        "warn_hits",
        "positive_hits",
        "negative_hits",
        "sign_consistency",
        "dataset_sign_consistency",
        "mean_effect",
        "median_effect",
        "mean_abs_effect",
        "max_abs_effect",
        "conserved_response_score",
        "datasets",
    ]
    summary = _write_csv(
        summary_rows,
        output_dir / "meta-summary.csv",
        ["rank", *gene_fieldnames],
    )
    consistency = _write_csv(
        gene_rows,
        output_dir / "gene-consistency.csv",
        gene_fieldnames,
    )
    dataset_summary = _write_csv(
        dataset_rows,
        output_dir / "dataset-meta-summary.csv",
        [
            "dataset_id",
            "input_contrasts",
            "included_contrasts",
            "input_top_gene_hits",
            "included_top_gene_hits",
        ],
    )

    result = {
        "schema_version": "1.0",
        "release_id": manifest.get("release_id"),
        "effects_manifest": str(effects_manifest_path),
        "effects_model": manifest.get("model"),
        "included_statuses": sorted(statuses),
        "top_n": top_n,
        "input_top_gene_hits": total_input_hits,
        "included_top_gene_hits": total_used_hits,
        "gene_count": len(gene_rows),
        "summary": summary,
        "gene_consistency": consistency,
        "dataset_summary": dataset_summary,
        "top_genes": summary_rows[:25],
        "datasets": dataset_rows,
    }
    _write_json(result, output_dir / "meta-manifest.json")
    return result
