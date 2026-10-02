# Metadata conventions

Keep four identifiers distinct: `source` (retrieval/annotation provider),
`study_id` (underlying study, e.g. a GEO accession), `sample_id` and
`patient_id`. The same study appearing in two portals remains one biological
cohort, not independent validation. Record portal-specific annotations and
their provenance separately.

Response must include its endpoint and source definition (for example RECIST
versus pathological complete response). Timepoints retain raw labels plus an
explicit ordered mapping. Cell type assignments include the published label
and a harmonized compartment; avoid claiming CNV-validated malignancy unless
that validation was actually performed. Study-specific values belong in the
run's metadata/manifest or an example mapping file.
