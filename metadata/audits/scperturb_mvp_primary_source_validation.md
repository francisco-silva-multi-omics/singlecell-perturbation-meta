# scPerturb MVP primary-source validation

Audit date: 2026-06-20

## Scope

This audit validates the five locally downloaded scPerturb RNA v1.4 H5AD files against their primary GEO/SRA records and the published scPerturb conversion code. It does not claim that a harmonized H5AD is byte-identical to an upstream file. Validation covers:

- study, sample, BioProject, and SRA accession lineage;
- local H5AD matrix dimensions, integer-count representation, and feature order;
- source-side matrix, barcode, gene, and identity-table dimensions where available;
- transformations made during scPerturb conversion;
- availability of controls and replicate/batch metadata;
- unresolved issues that must gate downstream analysis.

All five local `X` matrices contain integer-valued observations. None contains a separate `raw` object or counts layer, so the project must explicitly record that raw counts are stored in `X` for these releases.

## Overall result

| Local file | Primary record | Lineage result | Analysis admission |
|---|---|---|---|
| `AdamsonWeissman2016_GSM2406681_10X010.h5ad` | GSE90546 / GSM2406681 / SRX2400171 | Pass: exact source matrix shape, nonzero count, gene order, and transformed barcode order | Block pending perturbation/control curation |
| `DixitRegev2016_K562_TFs_7_days.h5ad` | GSE90063 / GSM2396858 / SRX2360555 | Pass: exact source matrix shape, nonzero count, gene order, and cell identities | Conditional: remove unassigned `nan` cells and repair `nperts` |
| `NormanWeissman2019_filtered.h5ad` | GSE133344 / PRJNA551220 / SRP212114 | Pass as a documented identity-annotated subset, not as the complete GEO filtered matrix | Pass with subset provenance retained |
| `SrivatsanTrapnell2020_sciplex3.h5ad` | GSE139944 / GSM4150378 / SRX7101188 | Pass: source pData/gene-defined dimensions and exact Ensembl order | Conditional: remove unassigned treatment hashes |
| `ReplogleWeissman2022_K562_essential.h5ad` | PRJNA831566 / SRP376262 / SAMN28561243 | Pass: primary processed H5AD shape, cell IDs, gene IDs, and every matrix value match | Pass with batch-aware analysis |

## Dataset findings

### Adamson et al. 2016

Primary evidence:

- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90546
- Exact sample: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2406681
- BioProject: PRJNA354963
- SRA study: SRP094654
- SRA experiment: SRX2400171

Reconciliation:

```text
GEO Matrix Market: 32,738 genes x 65,337 cells; 237,812,947 nonzero entries
Local H5AD:        65,337 cells x 32,738 genes; 237,812,947 nonzero entries
Gene IDs/order:    exact match
Cell order:        exact after the documented barcode suffix removal and make-unique step
```

The source has 65,257 cell-identity rows for 65,337 count-matrix barcodes. scPerturb removes the 10x gem-group suffix, creating 2,539 duplicate barcode roots, and deduplicates identity metadata before joining. The resulting H5AD has 2,613 cells with missing perturbation labels. Its apparent control guides are not normalized to `control`, so a control group cannot be selected safely from the harmonized field alone.

Decision: do not admit this file to perturbation-effect estimation until guide identities are rejoined using the original full barcode/gem-group identity and control guides are curated from the paper/source metadata.

### Dixit et al. 2016

Primary evidence:

- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90063
- Exact sample: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2396858
- BioProject: PRJNA354362
- SRA study: SRP093670
- SRA experiment: SRX2360555

Reconciliation:

```text
GEO Matrix Market: 23,111 genes x 33,013 cells; 121,279,430 nonzero entries
Local H5AD:        33,013 cells x 23,111 genes; 121,279,430 nonzero entries
Gene IDs/order:    exact match
Cell IDs/order:    exact match
```

The H5AD contains 5,381 cells normalized to `control`, based on intergenic guides. Another 4,979 cells have the literal perturbation value `nan`, while `nperts` is incorrectly set to 1 for every cell.

Decision: lineage passes. Exclude the 4,979 unassigned cells and recompute `nperts` before analysis.

### Norman et al. 2019

Primary evidence:

- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344
- BioProject: PRJNA551220
- SRA study: SRP212114

Reconciliation:

```text
GEO filtered matrix: 33,694 genes x 111,668 cells; 362,199,631 nonzero entries
GEO identity table:  111,445 cells
Local H5AD:          111,445 cells x 33,694 genes; 361,582,621 nonzero entries
Gene IDs/order:      exact match
Cell IDs/order:      exact identity-table match after barcode suffix removal and make-unique
```

The H5AD is an intentional subset of the GEO filtered matrix: it retains the 111,445 cells with identity annotations and drops 223 other filtered barcodes. It includes 11,855 normalized control cells. The original experiment is spread across multiple gem groups/libraries; their interpretation as biological versus technical replicates must follow the study design rather than treating cells as independent replicates.

Decision: pass as a documented subset. Record the 223-cell exclusion in derived manifests.

### Srivatsan et al. 2020, sci-Plex3

Primary evidence:

- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139944
- Exact sample: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4150378
- BioProject: PRJNA587707
- SRA study: SRP228591
- SRA experiment: SRX7101188

Reconciliation:

```text
GEO pData rows:       799,317 cells
GEO gene rows:        110,983 genes
Local H5AD:           799,317 cells x 110,983 genes
Ensembl IDs/order:    exact after documented version-suffix removal
```

The scPerturb conversion constructs the sparse UMI matrix directly from the GEO triplet count file using the pData and gene-table dimensions. The H5AD contains two named replicates and 17,578 vehicle/control cells. It also contains 36,522 cells with missing perturbation, replicate, plate, time, and dose values, consistent with unassigned treatment hashes.

Decision: lineage passes. Exclude the 36,522 unassigned cells from treatment comparisons; preserve replicate, plate, dose, and time during control matching.

### Replogle et al. 2022, K562 essential-scale

Primary evidence:

- Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/
- BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA831566
- SRA study: SRP376262
- K562 day-6 essential-scale BioSample: SAMN28561243
- Processed-data portal: https://gwps.wi.mit.edu/

The paper states that raw Perturb-seq sequencing is deposited under PRJNA831566 and processed single-cell data are distributed through the authors' portal. The SRA inventory contains 7,883 runs across four samples, including 624 runs assigned to `K562_day_6_essential_scale` / SAMN28561243.

The scPerturb conversion code reads `K562_essential_raw_singlecell_01.h5ad`, adds harmonized metadata, renames guide/batch fields, and calculates QC summaries. The local result has 310,385 cells, 8,563 genes, integer-valued `X`, 48 batch labels, and 10,691 controls.

The primary `K562_essential_raw_singlecell_01.h5ad` was obtained from the paper's Figshare deposit (file ID `35773219`). It is 10,661,879,995 bytes with file SHA-256 `3e5a63a9e892b21029bb55fca4e12517a49aad7af6c14133ca63d12cf68c6cee`.

Full row-streamed comparison results:

```text
Matrix shape:                 310,385 cells x 8,563 genes, exact
Cell identifiers/order:       exact
Ensembl identifiers/order:    exact
Every logical matrix value:   exact
Primary logical SHA-256:      1cb2ec238cec710df8ccfe612c3d4b661408b3b72e2e61ea931fb6f126ed7400
Harmonized logical SHA-256:   1cb2ec238cec710df8ccfe612c3d4b661408b3b72e2e61ea931fb6f126ed7400
```

Decision: full primary-source validation passes. The harmonized file changes metadata but preserves the primary count matrix, cells, and features exactly. Batch labels still require correct experimental-unit interpretation.

## Cross-dataset admission rules

1. Store `counts_location = X` explicitly for these five releases; do not infer it from the absence of a layer.
2. Exclude missing/unassigned perturbations before defining condition groups.
3. Curate control classes from primary metadata. A string heuristic is not sufficient for Adamson.
4. Treat cells as observations, not independent experimental replicates. Preserve gem group, batch, plate, and replicate fields and classify their experimental-unit meaning per paper.
5. Store primary accessions and this audit status in the dataset catalog.
6. Keep the scPerturb content hash/receipt as distribution provenance and the GEO/SRA accessions as experimental provenance; they serve different purposes.

## Conversion-code evidence

- scPerturb repository: https://github.com/sanderlab/scPerturb
- Adamson conversion: https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/scripts/AdamsonWeissman2016.py
- Dixit conversion: https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/scripts/DixitRegev2016.py
- Replogle conversion: https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/scripts/ReplogleWeissman2022.py
- sci-Plex conversion: https://github.com/sanderlab/scPerturb/blob/master/dataset_processing/scripts/SrivatsanTrapnell2020.py
