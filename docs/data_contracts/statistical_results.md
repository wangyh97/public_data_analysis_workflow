# Statistical result contract (v1)

One row is one explicitly defined hypothesis. Required: `analysis_id`,
`source`, `study_id`, `comparison`, `feature`, `method`, `n`, `estimate`,
`p_value`, `adjusted_p`, `adjustment_family`, `status`. Optional: confidence
interval, group sizes, paired count, effect unit, covariates and endpoint.

`n` is the number of independent observations used by the test. A paired test
also records `n_pairs`, pairing key and the exact baseline/follow-up labels.
`adjustment_family` states which hypotheses entered multiple-testing
correction. Plots display these stored values; they do not rerun tests.

If a test is impossible (for example insufficient pairs or constant values),
set the numerical fields to missing and give a machine-readable `status`.
