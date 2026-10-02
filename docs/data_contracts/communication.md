# Cell communication contract (v1)

One sample-level row contains `source`, `study_id`, `patient_id`, `sample_id`,
`timepoint`, `sender`, `receiver`, `pathway`, `strength`, `method`,
`database_version` and `n_cells`. Ligand-receptor edges also carry `ligand`,
`receptor` and an edge-level significance field. Patient-level summaries add
`n_samples` and state the within-patient aggregation rule.

Communication strength is model output, not a direct molecular measurement.
Changing database, downsampling, bootstraps or cell annotations changes its
meaning and requires a distinct processing version. When high/low groups are
formed from expression, record the split rule and patient assignment table.
