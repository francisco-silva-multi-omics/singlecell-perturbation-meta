from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import html
import json
import math
from pathlib import Path
from typing import Any

import h5py
import numpy as np


class PlotError(RuntimeError):
    """Raised when plot inputs are missing or inconsistent."""


COLORS = {
    "adamson": "#4C78A8",
    "dixit": "#F58518",
    "norman": "#54A24B",
    "replogle": "#B279A2",
    "sciplex3": "#E45756",
    "pass": "#54A24B",
    "warn": "#F2BE42",
    "fail": "#E45756",
    "control": "#4C78A8",
}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise PlotError(f"Expected JSON object: {path}")
    return payload


def _resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def _svg(width: int, height: int, body: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">\n'
        '<style>'
        'text{font-family:Arial,Helvetica,sans-serif;fill:#222}'
        '.title{font-size:22px;font-weight:700}'
        '.subtitle{font-size:12px;fill:#555}'
        '.axis{font-size:11px;fill:#555}'
        '.label{font-size:12px}'
        '.small{font-size:10px;fill:#555}'
        '.panel{fill:#fff;stroke:#ddd;stroke-width:1}'
        '</style>\n'
        '<rect width="100%" height="100%" fill="#ffffff"/>\n'
        f"{body}\n</svg>\n"
    )


def _text(x: float, y: float, value: Any, klass: str = "label", anchor: str = "start") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" class="{klass}" text-anchor="{anchor}">'
        f"{html.escape(str(value))}</text>"
    )


def _rect(x: float, y: float, width: float, height: float, fill: str, opacity: float = 1.0) -> str:
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(width, 0):.1f}" '
        f'height="{max(height, 0):.1f}" fill="{fill}" opacity="{opacity:.3f}"/>'
    )


def _line(x1: float, y1: float, x2: float, y2: float, stroke: str = "#999") -> str:
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="1"/>'


def _fmt_int(value: int | float) -> str:
    return f"{int(value):,}"


def _fmt_optional_int(value: int | float | None) -> str:
    if value is None:
        return "n/a"
    return _fmt_int(value)


def _bar_panel(
    *,
    x: int,
    y: int,
    width: int,
    height: int,
    title: str,
    values: list[tuple[str, float]],
    value_label: str,
) -> str:
    body = [_rect(x, y, width, height, "#fff"), _text(x + 14, y + 24, title, "label")]
    plot_x = x + 120
    plot_y = y + 44
    plot_w = width - 160
    row_h = min(28, (height - 64) / max(len(values), 1))
    max_value = max((value for _, value in values), default=1)
    for index, (name, value) in enumerate(values):
        yy = plot_y + index * row_h
        color = COLORS.get(name, "#777")
        body.append(_text(x + 14, yy + row_h * 0.65, name, "axis"))
        bar_w = 0 if max_value == 0 else (value / max_value) * plot_w
        body.append(_rect(plot_x, yy + 4, bar_w, row_h - 8, color))
        body.append(_text(plot_x + bar_w + 6, yy + row_h * 0.65, _fmt_int(value), "small"))
    body.append(_text(x + width - 14, y + height - 10, value_label, "small", "end"))
    return "\n".join(body)


def _load_dataset_summary(release_path: Path, effects_path: Path) -> list[dict[str, Any]]:
    release = _read_json(release_path)
    effects = _read_json(effects_path)
    by_dataset = {item["dataset_id"]: item for item in effects["datasets"]}
    rows = []
    for item in release["datasets"]:
        dataset_id = item["dataset_id"]
        effect = by_dataset.get(dataset_id, {})
        rows.append(
            {
                "dataset_id": dataset_id,
                "cells": int(item.get("cell_counts", {}).get("profile_n_obs") or 0),
                "genes": int(item.get("cell_counts", {}).get("profile_n_vars") or 0),
                "contrasts": int(effect.get("contrasts") or 0),
                "controls": int(effect["controls"]) if "controls" in effect else None,
                "status_counts": _load_effect_status_counts(effect),
            }
        )
    return rows


def _load_effect_status_counts(effect: dict[str, Any]) -> dict[str, int]:
    contrast_table = effect.get("contrast_table", {})
    path_value = contrast_table.get("path")
    if not path_value:
        return {}
    path = _resolve_path(path_value)
    if not path.exists():
        return {}
    counts: Counter[str] = Counter()
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "status" not in (reader.fieldnames or []):
            return {}
        for row in reader:
            status = row.get("status", "")
            if status:
                counts[status] += 1
    return dict(counts)


def _load_qc_counts(qc_summary_path: Path) -> tuple[dict[str, int], dict[str, int]]:
    qc = _read_json(qc_summary_path)
    return qc.get("dataset_status_counts", {}), qc.get("perturbation_status_counts", {})


def _plot_overview(
    *,
    rows: list[dict[str, Any]],
    qc_summary_path: Path,
    output_path: Path,
) -> None:
    _, perturbation_counts = _load_qc_counts(qc_summary_path)
    body = [
        _text(40, 38, "scPerturb MVP Pipeline Summary", "title"),
        _text(40, 58, "Cells, contrasts, and perturbation QC status after curation/harmonization/effect calculation.", "subtitle"),
    ]
    body.append(
        _bar_panel(
            x=40,
            y=80,
            width=520,
            height=230,
            title="Cells per dataset",
            values=[(row["dataset_id"], row["cells"]) for row in rows],
            value_label="cell count",
        )
    )
    body.append(
        _bar_panel(
            x=600,
            y=80,
            width=520,
            height=230,
            title="Effect contrasts per dataset",
            values=[(row["dataset_id"], row["contrasts"]) for row in rows],
            value_label="treated-vs-control contrasts",
        )
    )
    status_values = [(key, perturbation_counts.get(key, 0)) for key in ("pass", "warn", "fail", "control")]
    body.append(
        _bar_panel(
            x=40,
            y=340,
            width=520,
            height=210,
            title="Perturbation QC labels",
            values=status_values,
            value_label="labels",
        )
    )
    total_cells = sum(row["cells"] for row in rows)
    total_contrasts = sum(row["contrasts"] for row in rows)
    stats = [
        ("datasets", len(rows)),
        ("cells", _fmt_int(total_cells)),
        ("features max", _fmt_int(max(row["genes"] for row in rows))),
        ("effect contrasts", _fmt_int(total_contrasts)),
        ("QC pass labels", _fmt_int(perturbation_counts.get("pass", 0))),
        ("QC fail labels", _fmt_int(perturbation_counts.get("fail", 0))),
    ]
    body.append(_rect(600, 340, 520, 210, "#fff"))
    body.append(_text(614, 364, "MVP artifact counts", "label"))
    for index, (label, value) in enumerate(stats):
        yy = 398 + index * 24
        body.append(_text(630, yy, label, "axis"))
        body.append(_text(1020, yy, value, "label", "end"))
    _write(output_path, _svg(1160, 590, "\n".join(body)))


def _effect_stats(effects_manifest_path: Path) -> list[dict[str, Any]]:
    manifest = _read_json(effects_manifest_path)
    rows = []
    for dataset in manifest["datasets"]:
        matrix_path = _resolve_path(dataset["matrix"]["path"])
        with h5py.File(matrix_path, "r") as handle:
            delta = handle["delta_log1p_cpm"]
            contrast_count = int(delta.shape[0])
            feature_count = int(delta.shape[1])
            per_contrast_max_chunks = []
            per_contrast_p95_chunks = []
            for start in range(0, contrast_count, 128):
                chunk = np.abs(delta[start : start + 128])
                if chunk.size == 0:
                    continue
                per_contrast_max_chunks.append(chunk.max(axis=1))
                per_contrast_p95_chunks.append(np.percentile(chunk, 95, axis=1))
        if per_contrast_max_chunks:
            per_contrast_max = np.concatenate(per_contrast_max_chunks)
            per_contrast_p95 = np.concatenate(per_contrast_p95_chunks)
            median_max_abs = float(np.median(per_contrast_max))
            p90_max_abs = float(np.percentile(per_contrast_max, 90))
            median_p95_abs = float(np.median(per_contrast_p95))
        else:
            median_max_abs = 0.0
            p90_max_abs = 0.0
            median_p95_abs = 0.0
        rows.append(
            {
                "dataset_id": dataset["dataset_id"],
                "contrasts": contrast_count,
                "features": feature_count,
                "median_max_abs": median_max_abs,
                "p90_max_abs": p90_max_abs,
                "median_p95_abs": median_p95_abs,
            }
        )
    return rows


def _plot_effect_magnitude(rows: list[dict[str, Any]], output_path: Path) -> None:
    width, height = 980, 520
    margin_l, margin_b, margin_t, margin_r = 90, 80, 70, 40
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b
    max_y = max(row["p90_max_abs"] for row in rows) * 1.15
    body = [
        _text(40, 38, "Effect Magnitude by Dataset", "title"),
        _text(40, 58, "Per-contrast absolute log1p(CPM) delta summaries from treated-minus-control pseudobulk effects.", "subtitle"),
        _line(margin_l, margin_t, margin_l, margin_t + plot_h),
        _line(margin_l, margin_t + plot_h, margin_l + plot_w, margin_t + plot_h),
    ]
    for tick in range(6):
        value = max_y * tick / 5
        y = margin_t + plot_h - (value / max_y) * plot_h
        body.append(_line(margin_l - 5, y, margin_l + plot_w, y, "#e5e5e5"))
        body.append(_text(margin_l - 10, y + 4, f"{value:.2f}", "axis", "end"))
    bar_w = plot_w / (len(rows) * 3)
    for index, row in enumerate(rows):
        center = margin_l + (index + 0.5) * plot_w / len(rows)
        color = COLORS.get(row["dataset_id"], "#777")
        values = [
            ("median max", row["median_max_abs"], 0.95),
            ("p90 max", row["p90_max_abs"], 0.55),
            ("median p95", row["median_p95_abs"], 0.30),
        ]
        for offset, (_, value, opacity) in enumerate(values):
            x = center - bar_w * 1.5 + offset * bar_w
            h = (value / max_y) * plot_h
            body.append(_rect(x, margin_t + plot_h - h, bar_w * 0.85, h, color, opacity))
        body.append(_text(center, margin_t + plot_h + 22, row["dataset_id"], "axis", "middle"))
        body.append(_text(center, margin_t + plot_h + 38, f"n={row['contrasts']}", "small", "middle"))
    legend_x = width - 260
    for index, (label, opacity) in enumerate((("median max", 0.95), ("p90 max", 0.55), ("median p95", 0.30))):
        yy = 90 + index * 22
        body.append(_rect(legend_x, yy - 12, 16, 12, "#555", opacity))
        body.append(_text(legend_x + 24, yy, label, "small"))
    _write(output_path, _svg(width, height, "\n".join(body)))


def _top_gene_recurrence(effects_manifest_path: Path, *, top_n: int = 25) -> list[dict[str, Any]]:
    manifest = _read_json(effects_manifest_path)
    gene_scores: dict[str, dict[str, Any]] = defaultdict(lambda: {"datasets": set(), "hits": 0, "max_abs": 0.0})
    for dataset in manifest["datasets"]:
        dataset_id = dataset["dataset_id"]
        top_path = _resolve_path(dataset["top_genes"]["path"])
        with gzip.open(top_path, "rt", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                gene = row["gene_id"]
                effect = abs(float(row["effect_size"]))
                score = gene_scores[gene]
                score["datasets"].add(dataset_id)
                score["hits"] += 1
                score["max_abs"] = max(score["max_abs"], effect)
    rows = [
        {
            "gene_id": gene,
            "dataset_count": len(values["datasets"]),
            "hit_count": int(values["hits"]),
            "max_abs_effect": float(values["max_abs"]),
        }
        for gene, values in gene_scores.items()
    ]
    rows.sort(key=lambda row: (row["dataset_count"], row["hit_count"], row["max_abs_effect"]), reverse=True)
    return rows[:top_n]


def _plot_top_gene_recurrence(rows: list[dict[str, Any]], output_path: Path) -> None:
    width, height = 900, 760
    margin_l, margin_t = 180, 70
    plot_w = 630
    row_h = 24
    max_hits = max(row["hit_count"] for row in rows) if rows else 1
    body = [
        _text(40, 38, "Most Recurrent Top-Effect Genes", "title"),
        _text(40, 58, "Genes appearing most often among top absolute effects across all contrasts and datasets.", "subtitle"),
    ]
    for index, row in enumerate(rows):
        y = margin_t + index * row_h
        width_hits = (row["hit_count"] / max_hits) * plot_w
        color = "#4C78A8" if row["dataset_count"] >= 2 else "#9ecae9"
        body.append(_text(margin_l - 10, y + 16, row["gene_id"], "axis", "end"))
        body.append(_rect(margin_l, y + 4, width_hits, row_h - 7, color))
        body.append(_text(margin_l + width_hits + 6, y + 16, f"{row['hit_count']} hits, {row['dataset_count']} datasets", "small"))
    _write(output_path, _svg(width, height, "\n".join(body)))


def _write_report(
    *,
    output_path: Path,
    rows: list[dict[str, Any]],
    effect_rows: list[dict[str, Any]],
    recurrent_genes: list[dict[str, Any]],
    figures: dict[str, str],
) -> None:
    lines = [
        "# scPerturb MVP Visual Summary",
        "",
        "## Dataset Summary",
        "",
        "| Dataset | Cells | Genes | Contrasts | Pass | Warn | Controls |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        status_counts = row.get("status_counts", {})
        lines.append(
            f"| {row['dataset_id']} | {_fmt_int(row['cells'])} | {_fmt_int(row['genes'])} | "
            f"{_fmt_int(row['contrasts'])} | {_fmt_int(status_counts.get('pass', 0))} | "
            f"{_fmt_int(status_counts.get('warn', 0))} | {_fmt_optional_int(row['controls'])} |"
        )
    lines.extend(["", "## Effect Magnitude", "", "| Dataset | Contrasts | Median Max Abs Effect | P90 Max Abs Effect | Median P95 Abs Effect |", "|---|---:|---:|---:|---:|"])
    for row in effect_rows:
        lines.append(
            f"| {row['dataset_id']} | {row['contrasts']} | {row['median_max_abs']:.3f} | "
            f"{row['p90_max_abs']:.3f} | {row['median_p95_abs']:.3f} |"
        )
    lines.extend(["", "## Recurrent Top-Effect Genes", "", "| Gene | Datasets | Top-Gene Hits | Max Abs Effect |", "|---|---:|---:|---:|"])
    for row in recurrent_genes[:15]:
        lines.append(
            f"| {row['gene_id']} | {row['dataset_count']} | {row['hit_count']} | {row['max_abs_effect']:.3f} |"
        )
    lines.extend(["", "## Figures", ""])
    for label, path in figures.items():
        lines.append(f"- {label}: `{path}`")
    _write(output_path, "\n".join(lines) + "\n")


def plot_mvp_summary(
    *,
    release_manifest_path: Path,
    qc_summary_path: Path,
    effects_manifest_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    dataset_rows = _load_dataset_summary(release_manifest_path, effects_manifest_path)
    effect_rows = _effect_stats(effects_manifest_path)
    recurrent_genes = _top_gene_recurrence(effects_manifest_path)
    figures = {
        "overview": str(output_dir / "mvp_overview.svg"),
        "effect_magnitude": str(output_dir / "effect_magnitude.svg"),
        "top_gene_recurrence": str(output_dir / "top_gene_recurrence.svg"),
    }
    _plot_overview(rows=dataset_rows, qc_summary_path=qc_summary_path, output_path=output_dir / "mvp_overview.svg")
    _plot_effect_magnitude(effect_rows, output_dir / "effect_magnitude.svg")
    _plot_top_gene_recurrence(recurrent_genes, output_dir / "top_gene_recurrence.svg")
    _write_report(
        output_path=output_dir / "visual_summary.md",
        rows=dataset_rows,
        effect_rows=effect_rows,
        recurrent_genes=recurrent_genes,
        figures=figures,
    )
    manifest = {
        "schema_version": "1.0",
        "release_manifest": str(release_manifest_path),
        "qc_summary": str(qc_summary_path),
        "effects_manifest": str(effects_manifest_path),
        "effects_model": _read_json(effects_manifest_path).get("model", ""),
        "figures": figures,
        "report": str(output_dir / "visual_summary.md"),
        "datasets": dataset_rows,
        "effect_summary": effect_rows,
        "top_recurrent_genes": recurrent_genes,
    }
    _write(output_dir / "plot-manifest.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest
