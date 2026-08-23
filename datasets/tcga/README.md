# datasets/tcga -- TCGA dataset adapter

Implements the standard adapter interface for the **TCGA** public dataset
(UCSC Xena / Toil harmonized pan-cancer RNA-seq).

## Data source (chosen deliberately)

| Option | Verdict |
|---|---|
| **UCSC Xena Toil bulk** (used here) | One fixed-URL bulk file (~1.3 GB), checksummable, includes GTEx normals, no Bioconductor dependency, ideal for "download once, reuse forever". |
| GDC API | Authoritative and clinically richest, but paginated API, slow full-matrix pulls; reserved for a future `gdc_clinical` add-on. |
| TCGAbiolinks | Heavy Bioconductor stack, version-sensitive, brittle on HPC module R. Not used. |

Known limitation: Toil phenotype has limited clinical depth (no age/stage);
these standard fields are emitted as `NA` rather than fabricated. Add a GDC
clinical downloader later to enrich.

## Interface

| Script | Input | Output |
|---|---|---|
| `download.R --config config/datasets/tcga.yaml [--force]` | dataset YAML | `raw/TCGA/<files>` + `manifest.yaml` + `logs/download.log` |
| `prepare.R --config config/datasets/tcga.yaml [--cohorts SKCM,COAD] [--force]` | raw files | `processed/TCGA/<cohort>/{expression.rds, sample_metadata.rds/.tsv, feature_metadata.rds/.tsv, clinical.rds/.tsv, manifest.yaml}` + `logs/prepare.log` |
| `validate.R --processed-dir ... \| (--config ... --cohort ...)` | processed dir | PASS/FAIL report, exit 0/1 |

## Processed contract

* `expression.rds`: numeric matrix **feature x sample**, values `log2(TPM+0.001)` (passthrough).
* feature labels: gene symbols (Ensembl fallback flagged in `feature_metadata`);
  duplicate-symbol rows collapsed by configured policy (`max_iqr` default) with audit columns.
* one sample per patient by configurable priority (default Primary Tumor first).
* `clinical.rds`: standardized fields from the project-wide list; unavailable fields = `NA`
  (availability recorded in manifest).
* `manifest.yaml.fingerprint` = sha256 over preprocessing parameters + input checksums;
  unchanged fingerprint => prepare reuses existing outputs.

## Environment

```
export PDA_DATA_ROOT=/data/public_data
module load R/4.2.0-container        # HPC only; never inside these scripts
Rscript -e 'install.packages(c("yaml","data.table"))'   # data.table optional but strongly recommended
```
