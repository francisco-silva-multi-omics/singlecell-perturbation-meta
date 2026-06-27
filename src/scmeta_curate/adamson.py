from __future__ import annotations

import csv
import gzip
from pathlib import Path
from urllib.request import Request, urlopen

import h5py
import numpy as np

from scmeta_curate.h5ad import CurationError, dataframe_index


GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2406nnn/GSM2406681/suppl"
BARCODES_FILE = "GSM2406681_10X010_barcodes.tsv.gz"
IDENTITIES_FILE = "GSM2406681_10X010_cell_identities.csv.gz"
NEGATIVE_CONTROLS = {"63(mod)_pBA580", "Gal4-4(mod)_pBA582"}


def _fetch(url: str, path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    request = Request(url, headers={"User-Agent": "singlecell-perturbation-meta/0.1"})
    with urlopen(request, timeout=120) as response, temporary.open("wb") as handle:
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)
    temporary.replace(path)


def load_repair_metadata(input_path: Path, cache_dir: Path):
    barcode_path = cache_dir / BARCODES_FILE
    identities_path = cache_dir / IDENTITIES_FILE
    _fetch(f"{GEO_BASE}/{BARCODES_FILE}", barcode_path)
    _fetch(f"{GEO_BASE}/{IDENTITIES_FILE}", identities_path)

    with gzip.open(barcode_path, "rt", encoding="utf-8") as handle:
        source_barcodes = [line.strip() for line in handle if line.strip()]
    with gzip.open(identities_path, "rt", encoding="utf-8", newline="") as handle:
        identities = {row["cell BC"]: row for row in csv.DictReader(handle)}

    with h5py.File(input_path, "r") as h5ad:
        local_barcodes = dataframe_index(h5ad["obs"])
    if len(source_barcodes) != len(local_barcodes):
        raise CurationError("Adamson GEO barcode count does not match local H5AD")

    guides = np.empty(len(source_barcodes), dtype=object)
    reads = np.zeros(len(source_barcodes), dtype=np.int64)
    umis = np.zeros(len(source_barcodes), dtype=np.int64)
    controls = np.empty(len(source_barcodes), dtype=object)
    keep = np.zeros(len(source_barcodes), dtype=bool)
    for index, barcode in enumerate(source_barcodes):
        row = identities.get(barcode)
        if row is None:
            guides[index] = None
            controls[index] = None
            continue
        guide = row["guide identity"]
        keep[index] = True
        guides[index] = "control" if guide in NEGATIVE_CONTROLS else guide
        controls[index] = "negative_control" if guide in NEGATIVE_CONTROLS else "targeting"
        reads[index] = int(float(row["read count"]))
        umis[index] = int(float(row["UMI count"]))

    updates = {
        "source_cell_barcode": np.asarray(source_barcodes, dtype=object),
        "guide_id": guides.copy(),
        "perturbation": guides,
        "control_status": controls,
        "read count": reads,
        "UMI count": umis,
        "nperts": np.asarray(
            [0 if value == "negative_control" else 1 for value in controls], dtype=np.int64
        ),
    }
    return keep, updates
