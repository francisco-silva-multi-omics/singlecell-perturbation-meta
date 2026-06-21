# Single-cell perturbation meta-analysis

Reproducible acquisition and harmonization pipeline for a cross-study meta-analysis of single-cell genetic and pharmacological perturbations.

The initial implementation resolves immutable Zenodo releases into local manifests and downloads only explicitly selected artifacts. Downloads are resumable, checksum-verified, atomically promoted, and recorded in a receipt with a local SHA-256 digest.

## First data fetch

Resolve the scPerturb RNA v1.4 Zenodo record without downloading data:

```powershell
python scripts/fetch_data.py catalog zenodo `
  --record-id 13350497 `
  --output metadata/manifests/scperturb-rna-v1.4.json
```

Review the five-file MVP selection and transfer size:

```powershell
python scripts/fetch_data.py fetch `
  --manifest metadata/manifests/scperturb-rna-v1.4.json `
  --selection-file config/scperturb_mvp_files.json `
  --output-dir data/landing/scperturb/rna/v1.4 `
  --dry-run
```

Start the verified downloads after reviewing the dry run:

```powershell
python scripts/fetch_data.py fetch `
  --manifest metadata/manifests/scperturb-rna-v1.4.json `
  --selection-file config/scperturb_mvp_files.json `
  --output-dir data/landing/scperturb/rna/v1.4 `
  --workers 2 `
  --yes
```

Rerunning the command resumes `.partial` files and verifies already completed files. The default receipt is written beside the downloads as `download-receipt.json`.

To install the package and use the `scmeta-fetch` entry point:

```powershell
python -m pip install -e .
scmeta-fetch --help
```

## Tests

The test suite uses only the Python standard library and performs no network downloads:

```powershell
$env:PYTHONPATH = "src"
python -m unittest discover -s tests -v
```

See [DATABASE_AND_INGESTION_PLAN.md](DATABASE_AND_INGESTION_PLAN.md) for the catalog, storage, provenance, and scaling design.

