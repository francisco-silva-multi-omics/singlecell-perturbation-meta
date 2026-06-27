# Pipeline suite plan

This document defines the executable pipeline suite for the single-cell perturbation meta-analysis. It builds on the current fetch and curation tools and turns the project plan into staged, testable command groups.

## Design goals

1. Make every data transition reproducible from an immutable input artifact, a reviewed config, and a command-line invocation.
2. Keep raw public artifacts and landing H5AD files unchanged.
3. Treat metadata curation, control matching, QC, and perturbation-effect estimation as first-class pipeline stages with manifests.
4. Avoid one giant merged AnnData object. Use per-dataset curated matrices plus queryable metadata, summaries, pseudobulk tables, and release manifests.
5. Separate local execution from orchestration. Python CLIs perform the work; Snakemake can orchestrate the same commands once the stage contracts are stable.

## Current baseline

Implemented or already present in the working tree:

- `scmeta-fetch`: resolves Zenodo records and downloads selected files with checksums and receipts.
- `scmeta-curate`: filters/repairs the current scPerturb MVP cohort and writes a curation manifest.
- `metadata/audits/scperturb_mvp_primary_source_validation.md`: primary-source validation for the five MVP datasets.
- `data/landing/scperturb/rna/v1.4/download-receipt.json`: local acquisition receipt for the five H5AD files.
- `data/processed/scperturb/rna/v1.4/curation-manifest.json`: local curation manifest for Adamson, Dixit, sci-Plex3, Norman, and Replogle.

Known immediate gap:

- The Zenodo content-download fix and curation suite files are local changes and should be committed before larger pipeline work continues.

## Suite architecture

The suite should expose one top-level command namespace over focused packages:

```text
scmeta
  catalog       remote source manifests
  fetch         verified acquisition into landing storage
  audit         primary-source lineage validation
  curate        dataset-specific metadata repair and filtering
  release       define final analysis inputs and lineage
  profile       H5AD schema, counts, obs/var summaries, controls, batches
  harmonize     canonical obs/var/intervention tables
  qc            dataset and perturbation QC gates
  aggregate     replicate-aware pseudobulk and summary matrices
  contrast      explicit treated/control condition matching
  effects       perturbation effect estimation within study/context
  meta          cross-study meta-analysis
  grn           regulatory-network inference and TF-target validation
  drug          disease-reversal and pharmacological prioritization
  split         benchmark train/test split generation
  model         predictive model training and evaluation
  uncertainty   prediction intervals, calibration, and coverage
  report        HTML/Markdown reports and figure-ready tables
  db            optional PostgreSQL catalog migrations and loads
```

The existing `scmeta-fetch` and `scmeta-curate` entry points can remain as compatibility wrappers while the internal packages move toward this single namespace.

## Data layout

```text
config/
  sources/                 reviewed source intake records
  selections/              explicit file and dataset selections
  curation/                dataset-specific curation rules
  qc/                      QC thresholds and admission gates
  contrasts/               control matching and comparison definitions
  splits/                  benchmark split definitions

metadata/
  manifests/               remote source manifests
  audits/                  primary-source validation notes and JSON results
  profiles/                machine-readable H5AD profiles
  releases/                promoted release manifests
  db/                      SQL migrations and load receipts

data/
  raw/                     auxiliary primary files fetched from GEO, Figshare, etc.
  landing/                 immutable downloaded public artifacts
  processed/               curated per-dataset matrices
  harmonized/              canonical obs/var/intervention parquet and zarr stores
  qc/                      QC tables and reports
  aggregate/               pseudobulk and summary matrices
  results/                 effects, meta-analysis, model benchmarks
  networks/                regulatory networks and TF-target validation
  models/                  trained predictors and benchmark artifacts
  figures/                 generated manuscript figures
  reports/                 rendered analysis reports
```

Rule: every directory that contains derived artifacts needs a manifest or receipt that records inputs, code version, parameters, output sizes, hashes, and status.

## Pipeline stages

| Stage | Purpose | Primary inputs | Primary outputs | Status |
|---|---|---|---|---|
| 1. Intake | Curator-reviewed source registry | source YAML/JSON | normalized source records | planned |
| 2. Catalog | Resolve remote release metadata | intake records, Zenodo/GEO/Figshare APIs | source manifests | implemented for Zenodo |
| 3. Fetch | Download immutable artifacts | source manifests, selection files | landing files, receipt | implemented for Zenodo H5AD |
| 4. Audit | Validate lineage against primary records | landing files, primary support files | audit reports, admission decisions | partial/manual |
| 5. Curate | Repair metadata and remove invalid cells | landing files, curation rules | processed H5AD, curation manifest | partial for MVP |
| 6. Release | Define final analysis inputs and lineage | curation manifest, profiles, receipts, audits | release manifest | implemented for MVP |
| 7. Profile | Extract standardized H5AD metadata | release manifest | profiles, schema findings | implemented for MVP |
| 8. Harmonize | Create canonical analysis tables | release manifest, profiles, rules | obs/var/intervention/condition CSV | implemented for MVP |
| 9. QC | Gate datasets and perturbations | release manifest, profiles, harmonized metadata | QC reports, pass/fail tables | implemented for MVP |
| 10. Aggregate | Build pseudobulk summaries | harmonized metadata, final H5AD matrices | sparse pseudobulk matrices, group tables | implemented for intervention-level MVP |
| 11. Contrast | Match treated/control conditions | harmonized conditions | contrast definitions and admission status | implemented for condition and pooled MVP |
| 12. Effects | Estimate within-study perturbation effects | pseudobulk, contrasts | effect matrices and contrast tables | implemented for intervention and pooled MVP |
| 13. Meta | Combine comparable effects across studies | effect tables, study metadata | conserved/context-specific signatures | implemented for pooled MVP top-effect summaries |
| 14. GRN | Infer and validate regulatory networks | expression effects, optional ATAC, priors | TF-target networks, validation metrics | planned |
| 15. Drug | Score disease reversal and drug response | disease/drug signatures, dose/time metadata | drug rankings, reversal scores | planned |
| 16. Splits | Generate leakage-resistant train/test splits | harmonized metadata, effects | split manifests | planned |
| 17. Models | Train and score predictive models | matrices, splits, metadata | benchmark metrics, predictions | planned |
| 18. Uncertainty | Calibrate prediction confidence | predictions, bootstrap units | intervals, calibration tables | planned |
| 19. Reports | Produce reviewable outputs | manifests, QC, results | Markdown/HTML reports, figure tables | partial pooled MVP summary |
| 20. DB | Load catalog and summary metadata | manifests, profiles, QC, results | PostgreSQL control plane | planned |

## Stage contracts

### Intake

Input: reviewed records under `config/sources/*.json` or `*.yaml`.

Required fields:

- source name, version, DOI/accessions, release date, license if known;
- organism, modality, assay, perturbation type, biological context;
- expected files and expected count representation;
- primary-source links for experimental provenance;
- notes on controls, batches, replicates, and known ambiguities.

Output: normalized intake table and validation errors. No data files are downloaded in this stage.

### Catalog

Input: intake record plus remote provider.

Output: immutable manifest under `metadata/manifests/`.

For each artifact, store:

- provider, record ID, file name, URL, size, checksum, checksum type;
- source version and retrieval timestamp;
- selection tags such as `mvp`, `rna`, `genetic`, `drug`, `primary_support`.

Existing implementation covers Zenodo record `13350497` for scPerturb RNA v1.4.

### Fetch

Input: source manifest, selection file, output directory.

Output: landing files plus `download-receipt.json`.

Rules:

- use resumable `.partial` files;
- verify provider checksums where available;
- compute local SHA-256 for every artifact;
- never modify completed landing artifacts in place;
- fail if the selected byte count exceeds a configurable local budget unless explicitly approved.

### Audit

Input: landing H5AD files and primary support files from GEO/SRA/Figshare/project portals.

Output:

- human-readable audit report under `metadata/audits/`;
- optional JSON audit result per dataset;
- admission decision: `pass`, `conditional`, or `block`.

Minimum checks:

- matrix shape, nonzero count, and gene/cell order against primary source where available;
- counts representation: `X`, `raw`, or named layer;
- perturbation, control, batch, replicate, dose, and time availability;
- dataset-specific transformations made by upstream harmonization projects;
- unresolved confounding risks.

### Curate

Input: landing H5AD files, audit decisions, curation rules, auxiliary primary metadata.

Output: processed H5AD files and `curation-manifest.json`.

Rules:

- landing files remain immutable;
- all excluded cells are counted by reason;
- repaired columns preserve source fields when practical;
- control labels are curated from primary metadata, not inferred from loose string matching;
- passthrough datasets are still represented in the release manifest.

Current MVP behavior:

- Adamson: repairs guide/control metadata from GEO and removes cells without source identities.
- Dixit: removes literal `nan` perturbations and recomputes `nperts`.
- sci-Plex3: removes unassigned treatment hashes.
- Norman and Replogle: pass through with documented provenance.

### Release

Input: curation manifest, profile index, download receipt, and audit artifacts.

Output: release manifest under `metadata/releases/`.

The release manifest is the downstream contract. It resolves exactly one final analysis input per dataset and records:

- final H5AD path, size, SHA-256, and storage format;
- whether the file is curated or passthrough;
- curation counts and profile dimensions;
- lineage to landing input, download receipt, profile JSON, and audit files;
- admission status for analysis.

Current MVP release: `metadata/releases/scperturb_mvp_v0.1.json`.

### Profile

Input: release manifest.

Output: one JSON profile per dataset under `metadata/profiles/`.

Profile contents:

- matrix shape, encoding, dtype, sparsity, integer-count check;
- counts location and raw/layer availability;
- `obs` and `var` column inventories with encoding and missingness;
- perturbation/control counts;
- batch, replicate, dose, and time columns;
- high-level cell and gene QC distributions;
- warnings for ambiguous or missing required fields.

This is the next best stage to implement because it gives us machine-readable evidence for downstream harmonization and QC.

### Harmonize

Input: release manifest and final analysis H5AD files.

Output:

- canonical `obs` CSV per dataset;
- canonical `var` CSV per dataset;
- intervention table;
- condition table;
- optional AnnData Zarr per dataset for scalable reads.

Canonical cell-level fields:

```text
dataset_id
cell_id
organism
modality
assay
cell_context
cell_line
sample_id
experimental_unit_id
batch_id
replicate_id
perturbation_id
perturbation_type
perturbation_target
guide_id
drug_id
dose
dose_unit
time
time_unit
is_control
control_type
counts_location
source_cell_barcode
```

Do not force all datasets to have all fields. Missing values must be explicit and explained by profile warnings or audit notes.

Current MVP output: `data/harmonized/scperturb/rna/v0.1/`.

### QC

Input: profiles, harmonized metadata, processed matrices.

Output: dataset-level and perturbation-level QC tables.

Dataset gates:

- count matrix is recoverable and integer-like;
- perturbation and control assignments are unambiguous;
- at least one valid control group exists for intended contrasts;
- batch/replicate metadata is available or absence is documented;
- confounding between perturbation and batch/time/plate is quantified.

Perturbation gates:

- minimum cells per condition;
- minimum experimental units where available;
- control match exists;
- not fully aliased with a single technical batch unless the analysis model explicitly permits it.

Current MVP output: `metadata/qc/scperturb_mvp_v0.1/`.

### Aggregate

Input: processed matrices, harmonized metadata, contrast/control-match definitions.

Output:

- pseudobulk expression matrices by dataset, context, perturbation, and experimental unit;
- per-condition summary statistics;
- replicate and cell-count tables.

Aggregation must preserve experimental units. Cells are observations, not biological replicates.

Current MVP outputs:

- `data/aggregate/scperturb/rna/v0.1/intervention/`;
- `data/aggregate/scperturb/rna/v0.1/condition/`.

The first implemented grouping is `intervention`, which aggregates by perturbation identity and keeps the output tractable across all five datasets. The implementation also supports `bio-condition` and `condition`; those should be used selectively when dose/time or batch/replicate-aware analysis is needed.

### Contrast

Input: harmonized `conditions.csv`.

Output:

- `condition-contrasts.csv`;
- `contrast-summary.json`;
- `pooled-contrasts.csv`;
- `pooled-contrast-summary.json`.

Current MVP outputs:

- `data/results/scperturb/rna/v0.1/contrasts/condition/`;
- `data/results/scperturb/rna/v0.1/contrasts/pooled/`.

Controls are matched by dataset, cell line, cell type, disease, time, batch, and replicate when available. Dose is retained as treated metadata but excluded from control matching so drug-treated conditions can match vehicle controls at the same context/time/batch/replicate.

The pooled contrast mode groups treated conditions by dataset, biological context, perturbation, dose, and time, then requires matched controls in the same batch/replicate strata. This solves the sparse per-condition problem in CRISPR-scale screens while preserving replicate labels for downstream modeling.

### Effects

Input: pseudobulk matrices and contrast definitions.

Output: within-study effect matrices and contrast tables.

Current MVP outputs:

- `data/results/scperturb/rna/v0.1/effects/intervention/`;
- `data/results/scperturb/rna/v0.1/effects/pooled/`.

The pooled effect mode consumes condition-level pseudobulk and pooled contrast definitions. It sums matched treated and control condition rows per contrast, computes `log1p(CPM)` treated-minus-control deltas, and preserves contrast status plus stratum counts in the output contrast table.

Minimum fields:

```text
dataset_id
context_id
perturbation_id
target_id
gene_id
effect_size
log_fold_change
p_value
adjusted_p_value
n_cells_treated
n_cells_control
n_units_treated
n_units_control
model
covariates
qc_status
```

The first implementation can use a conservative pseudobulk differential-expression workflow before adding virtual-cell model benchmarks.

Current MVP output: `data/results/scperturb/rna/v0.1/effects/intervention/`.

The first implemented effect model is `treated_minus_control_delta` on intervention-level pseudobulk counts:

```text
effect_size = log1p(CPM_treated) - log1p(CPM_control)
```

This is an effect-size layer, not a statistical DE test. P-values, FDR, replicate-aware models, and covariates should be added after batch/replicate-aware pseudobulk is produced.

### Meta-analysis

Input: effect tables and comparable context/perturbation mappings.

Output:

- conserved perturbation signatures;
- context-specific signatures;
- heterogeneity statistics;
- perturbation-by-gene-by-context tables.

Only combine effects when the intervention, biological context, modality, and control definition are comparable enough to defend.

Current MVP output: `metadata/meta/scperturb_mvp_v0.1_pooled/`.

Implemented artifacts:

- `meta-summary.csv`;
- `gene-consistency.csv`;
- `dataset-meta-summary.csv`;
- `meta-manifest.json`.

The MVP meta stage summarizes pooled top-effect rows rather than fitting a full statistical meta-analysis model. It ranks genes by dataset coverage, hit count, sign consistency, and effect magnitude while keeping `pass` and `warn` evidence counts separate.

### Regulatory Networks

Input: effect tables, RNA expression matrices, optional ATAC/multiome features, and external priors.

Output:

- RNA-only regulatory networks;
- RNA+ATAC regulatory networks when multiome data is added;
- TF-target validation tables using observed CRISPRi/a perturbations;
- consensus network edges with evidence categories and weights.

First implementation:

- start RNA-only;
- use TF perturbations from Dixit/Norman/Replogle where applicable;
- test whether predicted TF targets move in the expected direction after perturbation;
- compare against simple correlation networks before adding complex GRN tools.

Later implementation:

- add SCENIC+/CellOracle/LINGER-style multiome networks;
- connect enhancer/promoter evidence when ATAC is available;
- add motif and pathway priors as explicit evidence columns.

### Drug Reversal

Input: disease signatures, drug signatures, harmonized dose/time metadata, and effect tables.

Output:

- reversal scores by cell context, disease state, drug, dose, and time;
- ranked drug and target tables;
- sensitive/resistant cell-context summaries;
- uncertainty-aware prioritization tables.

The score must be stratified. A single global disease/drug correlation is too coarse for this project.

Minimum formula:

```text
disease_signature = disease - healthy
drug_signature = drug_treated - untreated
reversal_score = -correlation(disease_signature, drug_signature)
```

Required output keys:

```text
cell_context
disease_status
drug_id
drug_name
dose
time
reversal_score
n_cells_treated
n_cells_control
n_units_treated
n_units_control
uncertainty_score
qc_status
```

### Splits

Input: harmonized metadata, processed matrices, effect tables.

Output: split definitions for downstream benchmarks.

Required split families:

- seen perturbation, seen context;
- unseen perturbation;
- unseen biological context;
- unseen dataset/study;
- combination generalization where available.

Leakage controls:

- no same-cell or same-condition leakage across train/test;
- split by perturbation and study when testing generalization;
- record all excluded perturbations and why.

### Models

Input: processed matrices, harmonized metadata, effect tables, and split manifests.

Output:

- model predictions;
- benchmark metrics;
- trained model artifacts where appropriate;
- model cards describing training data, split type, inputs, and limitations.

Baseline models are required before complex models:

```text
mean baseline
linear/ridge/elastic net
kNN
random forest or XGBoost
simple scVI latent model
```

Specialized models can then be evaluated:

```text
CPA
chemCPA
GEARS
scGen
scGPT or Geneformer fine-tuning
compact Transformer or SSM/Mamba hypothesis model
```

Do not treat complex models as useful unless they beat simple baselines under strict splits.

Metrics:

```text
Pearson/Spearman on delta-expression
RMSE on logFC
DE gene overlap
precision@k and recall@k
sign accuracy
pathway overlap
TF-target recovery
calibration error
uncertainty coverage
```

### Uncertainty

Input: benchmark predictions, effect tables, replicate/bootstrap units, and model outputs.

Output:

- prediction intervals by gene/pathway/drug;
- calibration curves or tables;
- uncertainty coverage metrics;
- uncertainty-aware rankings.

Methods to support progressively:

```text
bootstrap by study or experimental unit
deep ensembles
Monte Carlo dropout
Bayesian last layer
conformal prediction
```

### Reports

Input: receipts, profiles, QC tables, effects, meta-analysis outputs.

Output:

- dataset acquisition report;
- curation and admission report;
- QC report;
- perturbation-effect report;
- regulatory-network report;
- predictive-model benchmark report;
- drug-reversal prioritization report;
- figure-ready tables for the manuscript.

Reports should be generated from machine-readable artifacts, not manually copied notebook output.

### Database

The database is a metadata control plane, not an expression store.

Initial PostgreSQL tables:

- `source`, `source_release`, `artifact`, `acquisition`;
- `dataset`, `assay`, `biosample`, `experimental_unit`;
- `intervention`, `intervention_component`, `condition`, `control_match`;
- `artifact_lineage`, `profile_summary`, `qc_result`;
- `effect_result`, `meta_result`.

Expression matrices remain in H5AD/Zarr; cell-level metadata can be Parquet. PostgreSQL stores summaries, keys, provenance, and artifact pointers.

## Orchestration plan

Short term: keep direct CLI commands documented in README while stage contracts are still changing.

Medium term: add `workflow/Snakefile` with rules that call the CLIs:

```text
catalog -> fetch -> audit -> curate -> release -> profile -> harmonize -> qc -> aggregate -> effects -> meta -> grn -> drug -> split -> model -> uncertainty -> report
```

Each Snakemake rule should only depend on declared file inputs/outputs and config files. The Python CLIs must remain runnable outside Snakemake for debugging and tests.

## MVP implementation order

1. Commit current fetch fix and curation work.
2. Add `scmeta-profile` to inspect the five local H5AD files and write JSON profiles. Done for the MVP cohort.
3. Add `metadata/releases/scperturb_mvp_v0.1.json` that points to processed outputs and passthrough landing files. Done.
4. Add harmonized metadata export for the five datasets, starting with `obs` parquet or CSV if parquet dependencies are deferred. Done with compressed CSV for v0.1.
5. Add QC gates using profile and curation manifests.
6. Add contrast/control-match config for the first defensible comparisons.
7. Build pseudobulk aggregation for one genetic dataset and one drug dataset.
8. Add first perturbation-effect tables.
9. Add strict benchmark split manifests.
10. Add a baseline model benchmark before any foundation-model work.
11. Add RNA-only TF-target recovery analysis for TF perturbation datasets.
12. Add first drug-reversal scoring on sci-Plex3 or another drug dataset.
13. Add uncertainty estimates for effects and predictions.
14. Add PostgreSQL schema only after the metadata contracts above stabilize.
15. Add Snakemake orchestration after at least `fetch -> curate -> profile -> qc` are stable.

## Near-term engineering tasks

### Task 1: profile command

Create `src/scmeta_profile/` and `scripts/profile_data.py`.

Command sketch:

```powershell
python scripts/profile_data.py profile `
  --curation-manifest data/processed/scperturb/rna/v1.4/curation-manifest.json `
  --output-dir metadata/profiles/scperturb/rna/v1.4
```

The command produces one JSON per dataset plus an index file. The MVP profiles are written under `metadata/profiles/scperturb/rna/v1.4`.

### Task 2: release manifest

Create a release manifest that resolves the final analysis input for each MVP dataset:

- Adamson: processed H5AD;
- Dixit: processed H5AD;
- sci-Plex3: processed H5AD;
- Norman: passthrough landing H5AD;
- Replogle: passthrough landing H5AD.

This avoids pretending that passthrough datasets were physically copied into `data/processed`.

### Task 3: QC command

Create `src/scmeta_qc/`.

First checks:

- required fields present or documented missing;
- perturbation and control counts;
- missing perturbation count after curation;
- cell count by batch/replicate/condition;
- perturbation-batch contingency summaries.

### Task 4: first harmonized metadata export

Create canonical metadata tables without rewriting matrices.

Preferred output:

```text
data/harmonized/scperturb/rna/v0.1/obs/{dataset_id}.parquet
data/harmonized/scperturb/rna/v0.1/var/{dataset_id}.parquet
data/harmonized/scperturb/rna/v0.1/interventions.parquet
data/harmonized/scperturb/rna/v0.1/conditions.parquet
```

Fallback if we want to avoid extra dependencies temporarily: compressed CSV with stable schemas.

### Task 5: split manifest generator

Create strict benchmark split definitions before training any model.

Required split files:

```text
metadata/splits/seen_perturbation_seen_context.json
metadata/splits/unseen_perturbation.json
metadata/splits/unseen_context.json
metadata/splits/unseen_dataset.json
metadata/splits/unseen_combination.json
```

### Task 6: first effect table

Generate one conservative pseudobulk differential-expression result for a defensible genetic perturbation comparison.

The goal is not model performance yet. The goal is to prove that the metadata, controls, experimental units, and count matrix can support one end-to-end scientific contrast.

### Task 7: figure/report scaffolding

Create report outputs that map to the manuscript figures in the thesis plan:

```text
Figure 1: dataset and pipeline design
Figure 2: atlas of responses
Figure 3: conservation versus context
Figure 4: regulatory networks
Figure 5: model benchmark
Figure 6: causal regulatory recovery
Figure 7: pharmacological reversal
```

Current MVP pooled summary output: `figures/scperturb_mvp_v0.1_pooled/`.

Implemented artifacts:

- `mvp_overview.svg`;
- `effect_magnitude.svg`;
- `top_gene_recurrence.svg`;
- `visual_summary.md`;
- `plot-manifest.json`.

## Definition of done for suite v0.1

The suite v0.1 is complete when:

- the five MVP datasets can be fetched or verified from local receipts;
- curation can be rerun from landing files and primary support files;
- each final analysis artifact has a release manifest entry and SHA-256;
- each dataset has a JSON profile;
- each dataset and perturbation has a QC admission status;
- canonical metadata tables exist for cells, features, interventions, and conditions;
- at least one genetic and one drug perturbation contrast can be aggregated into pseudobulk;
- at least one within-study effect table is produced;
- strict split manifests exist before model training starts;
- at least one baseline model is evaluated;
- the commands are covered by unit tests and can be orchestrated from a single documented command path.

## Longer-term suite milestones

### v0.2: First Statistical Results

- QC-passing perturbation contrasts for the five MVP datasets.
- Pseudobulk differential expression for selected genetic and drug perturbations.
- Perturbation x gene x context effect table.
- First conservation/context-dependence report.

### v0.3: Benchmark Layer

- Strict train/test splits.
- Baseline model leaderboard.
- First specialized model comparison.
- Leakage audit for all benchmarks.

### v0.4: Regulatory Layer

- RNA-only GRN baseline.
- TF-target recovery validation using CRISPRi/a perturbations.
- Consensus regulatory edge table with evidence categories.

### v0.5: Pharmacological Layer

- Disease/drug signature schema.
- Stratified disease-reversal scores.
- Drug ranking with uncertainty.

### v1.0: Manuscript-Ready Pipeline

- Reproducible Snakemake workflow.
- Frozen data release manifest.
- Full QC, effects, meta-analysis, benchmark, GRN, drug, and uncertainty reports.
- Figure-ready tables for the manuscript deliverables described in the thesis plan.
