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

## Curate the MVP cohort

The curation command leaves landing files unchanged, repairs Adamson metadata from GEO, removes unassigned Dixit and sci-Plex3 cells, and writes checksummed outputs under `data/processed`:

```powershell
$env:PYTHONPATH = "src"
python scripts/curate_data.py curate `
  --input-dir data/landing/scperturb/rna/v1.4 `
  --output-dir data/processed/scperturb/rna/v1.4 `
  --geo-cache data/raw/geo/GSE90546/GSM2406681
```

Verify the primary Replogle processed matrix against the harmonized release:

```powershell
python scripts/curate_data.py verify-replogle `
  --primary data/raw/figshare/replogle/K562_essential_raw_singlecell_01.h5ad `
  --harmonized data/landing/scperturb/rna/v1.4/ReplogleWeissman2022_K562_essential.h5ad `
  --output metadata/audits/replogle_processed_verification.json
```

## Build the MVP release manifest

Create the manifest that defines the final analysis input for each MVP dataset:

```powershell
$env:PYTHONPATH = "src"
python scripts/release_data.py build `
  --release-id scperturb_mvp_v0.1 `
  --curation-manifest data/processed/scperturb/rna/v1.4/curation-manifest.json `
  --profile-index metadata/profiles/scperturb/rna/v1.4/profile-index.json `
  --download-receipt data/landing/scperturb/rna/v1.4/download-receipt.json `
  --primary-audit metadata/audits/scperturb_mvp_primary_source_validation.md `
  --replogle-verification metadata/audits/replogle_processed_verification.json `
  --output metadata/releases/scperturb_mvp_v0.1.json
```

The release manifest is the stable boundary for downstream QC and harmonization. It resolves curated outputs and passthrough landing files into one final H5AD per dataset, with lineage to receipts, profiles, and audits.

## Profile the curated release

Generate machine-readable H5AD profiles for the final MVP analysis inputs:

```powershell
$env:PYTHONPATH = "src"
python scripts/profile_data.py profile `
  --release-manifest metadata/releases/scperturb_mvp_v0.1.json `
  --output-dir metadata/profiles/scperturb/rna/v1.4
```

The command writes one `*.profile.json` per dataset and a `profile-index.json` summary. Profiles include matrix shape, encoding, count-location evidence, `obs`/`var` columns, perturbation/control/batch summaries, and warnings for missing analysis metadata.

## Run MVP QC gates

Run dataset-level and perturbation-level QC against the release manifest:

```powershell
$env:PYTHONPATH = "src"
python scripts/qc_data.py qc `
  --release-manifest metadata/releases/scperturb_mvp_v0.1.json `
  --output-dir metadata/qc/scperturb_mvp_v0.1 `
  --min-cells 50 `
  --min-control-cells 50
```

The QC stage writes `dataset-qc.json/csv`, `perturbation-qc.json/csv`, and `qc-summary.json`. The first gates check matrix recoverability, integer-like counts, perturbation labels, control availability, cell-count thresholds, and whether batch/replicate confounding can be assessed.

## Harmonize MVP metadata

Export canonical cell, feature, intervention, and condition metadata:

```powershell
$env:PYTHONPATH = "src"
python scripts/harmonize_data.py harmonize `
  --release-manifest metadata/releases/scperturb_mvp_v0.1.json `
  --output-dir data/harmonized/scperturb/rna/v0.1
```

The harmonization stage writes compressed per-dataset `obs` and `var` CSV files, plus release-level `interventions.csv`, `conditions.csv`, and `harmonization-manifest.json`.

## Aggregate MVP pseudobulk counts

Build intervention-level pseudobulk count matrices from the harmonized metadata and final H5AD matrices:

```powershell
$env:PYTHONPATH = "src"
python scripts/aggregate_data.py aggregate `
  --harmonization-manifest data/harmonized/scperturb/rna/v0.1/harmonization-manifest.json `
  --output-dir data/aggregate/scperturb/rna/v0.1 `
  --group-by intervention `
  --chunk-size 8192
```

The aggregation stage writes one sparse HDF5 pseudobulk matrix per dataset, a `groups.csv` file per dataset, and a `pseudobulk-manifest.json`. The default `intervention` grouping is the first tractable MVP layer; finer `bio-condition` and `condition` groupings are available for later dose/time and batch-aware analyses.

Build condition-level pseudobulk when running pooled contrast-aware effects:

```powershell
$env:PYTHONPATH = "src"
python scripts/aggregate_data.py aggregate `
  --harmonization-manifest data/harmonized/scperturb/rna/v0.1/harmonization-manifest.json `
  --output-dir data/aggregate/scperturb/rna/v0.1 `
  --group-by condition `
  --chunk-size 8192
```

## Compute MVP effect matrices

Compute treated-minus-control effect matrices from intervention-level pseudobulk counts:

```powershell
$env:PYTHONPATH = "src"
python scripts/effect_data.py effects `
  --pseudobulk-manifest data/aggregate/scperturb/rna/v0.1/intervention/pseudobulk-manifest.json `
  --output-dir data/results/scperturb/rna/v0.1/effects/intervention `
  --min-cells 50 `
  --top-genes 100
```

The first effect stage uses `log1p(CPM)` pseudobulk expression and writes treated-minus-control deltas. Outputs include one `effects.h5` matrix per dataset, one `contrasts.csv` per dataset, `top_genes.csv.gz` for inspection, and a release-level `effects-manifest.json`.

After building pooled contrasts and condition-level pseudobulk, compute matched pooled effects:

```powershell
$env:PYTHONPATH = "src"
python scripts/effect_data.py pooled `
  --pseudobulk-manifest data/aggregate/scperturb/rna/v0.1/condition/pseudobulk-manifest.json `
  --pooled-contrasts data/results/scperturb/rna/v0.1/contrasts/pooled/pooled-contrasts.csv `
  --output-dir data/results/scperturb/rna/v0.1/effects/pooled `
  --top-genes 100
```

The pooled effect stage includes `pass` and `warn` contrasts by default, excludes `fail`, and writes one effect matrix, contrast table, and top-gene table per dataset. In the MVP run it produced 4,613 effects with zero skipped contrasts.

## Summarize pooled meta-signatures

Rank recurrent top-effect genes by dataset coverage, sign consistency, hit count, and effect magnitude:

```powershell
$env:PYTHONPATH = "src"
python scripts/meta_data.py summarize `
  --effects-manifest data/results/scperturb/rna/v0.1/effects/pooled/effects-manifest.json `
  --output-dir metadata/meta/scperturb_mvp_v0.1_pooled `
  --top-n 100
```

The meta stage writes `meta-summary.csv`, `gene-consistency.csv`, `dataset-meta-summary.csv`, and `meta-manifest.json`. In the MVP run it summarized 24,359 genes from 461,300 top-effect hits.

## Build condition-level contrasts

Create explicit treated/control contrast definitions from harmonized conditions:

```powershell
$env:PYTHONPATH = "src"
python scripts/contrast_data.py build `
  --harmonization-manifest data/harmonized/scperturb/rna/v0.1/harmonization-manifest.json `
  --output-dir data/results/scperturb/rna/v0.1/contrasts/condition `
  --min-treated-cells 50 `
  --min-control-cells 50
```

The contrast stage writes `condition-contrasts.csv` and `contrast-summary.json`. Controls are matched by dataset, cell context, disease, time, batch, and replicate where available. Drug dose is retained on treated contrasts but not required for vehicle-control matching.

For sparse guide or perturbation strata, build pooled replicate-aware contrasts. This keeps batch/replicate as matched strata, but admits a perturbation/context contrast when the matched strata clear the pooled cell thresholds:

```powershell
$env:PYTHONPATH = "src"
python scripts/contrast_data.py build-pooled `
  --harmonization-manifest data/harmonized/scperturb/rna/v0.1/harmonization-manifest.json `
  --output-dir data/results/scperturb/rna/v0.1/contrasts/pooled `
  --min-treated-cells 50 `
  --min-control-cells 50 `
  --min-strata 2
```

The pooled stage writes `pooled-contrasts.csv` and `pooled-contrast-summary.json`. In the MVP run this produces 4,255 pass, 358 warn, and 262 fail pooled contrasts. Warnings indicate enough cells but only one matched stratum.

## Plot MVP summaries

Generate dependency-light SVG figures from release, QC, and pooled effect artifacts:

```powershell
$env:PYTHONPATH = "src"
python scripts/plot_results.py mvp-summary `
  --release-manifest metadata/releases/scperturb_mvp_v0.1.json `
  --qc-summary metadata/qc/scperturb_mvp_v0.1/qc-summary.json `
  --effects-manifest data/results/scperturb/rna/v0.1/effects/pooled/effects-manifest.json `
  --output-dir figures/scperturb_mvp_v0.1_pooled
```

The plotting stage writes `mvp_overview.svg`, `effect_magnitude.svg`, `top_gene_recurrence.svg`, `visual_summary.md`, and `plot-manifest.json`. The plotter streams large effect matrices in chunks, so the pooled sci-Plex3 matrix can be summarized without loading the full matrix into memory.
