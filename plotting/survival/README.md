# plotting/survival -- KM curves + Cox forest plot

Reads **only** the two tables produced by `analyses/survival`:

```text
INPUT
  --results results/<name>/tables/survival_results.tsv
  --curves  results/<name>/tables/km_curves.tsv
  [--plot-config config/plotting/default.yaml]
  [--out-dir results/<name>/figures] [--types km,forest] [--prefix survival]

OUTPUT
  survival_km.<pdf|png>       Kaplan-Meier step curves per block,
                              annotated with log-rank p/q read from the table
  survival_forest.<pdf|png>   HR (log scale) with 95% CI, faceted by model
```

No statistics are recomputed here; annotations are lookups from the results
table. Changing style only touches `config/plotting/default.yaml`.
