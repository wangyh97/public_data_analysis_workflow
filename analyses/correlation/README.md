# analyses/correlation -- correlation analysis module

Dataset-agnostic: reads **only** the standardized processed-data contract.
This module contains zero knowledge of TCGA, GDC, Xena or any API.

## Interface

```text
INPUT
  --config config/analyses/<analysis>.yaml
    analysis.name / analysis.type=correlation
    datasets: [{source, cohort}]            # resolved via dataset registry + PDA_DATA_ROOT
    targets.genes | targets.geneset_file
    comparison.genes | comparison.geneset_file
    correlation.method (spearman|pearson), adjust_method (BH|holm|none),
    adjust_family (target_gene|global), min_n

OUTPUT
  results/<analysis.name>/tables/correlation_results.tsv   long format:
      dataset, cohort, target_gene, comparison_gene, n,
      correlation, p_value, p_adjust, method, status, processed_fingerprint
  results/<analysis.name>/metadata/config_used.yaml        resolved input config echo
  results/<analysis.name>/metadata/analysis_metadata.yaml  provenance (inputs, software, timing)
  results/<analysis.name>/logs/run.log; sessionInfo.txt
```

## Behaviour notes

* Missing requested genes are **warnings**; the run fails only if none remain.
* `status = insufficient_n` rows keep NA statistics and are excluded from FDR families.
* p-adjustment is computed within `dataset x cohort x target_gene` by default.
* Self-pairs (gene in both lists) are skipped with a warning.
* Re-running always recomputes (cheap); plotting is a separate module.

## Example

```bash
export PDA_DATA_ROOT=/data/public_data
export PDA_RESULTS_ROOT=/project/results
Rscript analyses/correlation/run.R --config config/analyses/correlation_example.yaml
```
