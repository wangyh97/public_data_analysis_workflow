# analyses/survival -- survival analysis module

Dataset-agnostic: reads **only** the standardized processed contract
(`expression.rds`, `clinical.rds`, `manifest.yaml`). Zero knowledge of any
data source or API.

## Interface

```text
INPUT
  --config config/analyses/<analysis>.yaml
    datasets: [{source, cohort}]
    targets.genes | targets.geneset_file          # one KM/HR analysis per gene
    survival.endpoint        OS | PFS | DSS | DFI (columns in clinical.rds)
    survival.split_method    median | fixed (+ cutoff)
    survival.covariates      e.g. [age, sex]; absent fields warn and drop
    survival.models          logrank | cox_univariate | cox_adjusted
    survival.min_*           hard floors; below => loud failure (fail loudly)

OUTPUT
  results/<name>/tables/survival_results.tsv   long format:
      dataset, cohort, target_gene, endpoint, model, term,
      estimate_type (chi2|hazard_ratio), n, n_events,
      estimate, ci_low, ci_high, p_value, p_adjust, status
  results/<name>/tables/km_curves.tsv          sample/step data for plotting:
      dataset, cohort, target_gene, endpoint, group(High|Low),
      time, surv, n_risk, n_event, n_group, events_group,
      split_method, cutpoint
  results/<name>/metadata/{config_used,analysis_metadata}.yaml
  results/<name>/logs/run.log ; sessionInfo.txt
```

## Statistics notes

* Expression is z-scored within dataset x cohort; Cox HR is **per 1 SD**.
* High/Low split: `median` = strictly above median (deterministic);
  `fixed` = user cutoff on raw expression. No automatic "optimal cutpoint"
  by design (overfitting risk) -- add later only with explicit caveat.
* Cleaning rules (all logged): one sample per patient, time>0, event in {0,1};
  exclusions counted with reasons; below floors => error, never silent.
* FDR within dataset x cohort x target x model over tested rows.

## Example

```bash
Rscript analyses/survival/run.R --config config/analyses/survival_example.yaml
```
