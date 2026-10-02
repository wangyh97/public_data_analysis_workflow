# Patient observation contract (v1)

One row is a patient, treatment node, cell compartment and feature. Required:
`source`, `study_id`, `patient_id`, `timepoint`, `compartment`, `feature`,
`value`, `unit`, `n_cells`. The combination
`source + study_id + patient_id + timepoint + compartment + feature` is unique
within a processing version. Optional: `specimen_id`, `response`,
`response_endpoint`, `annotation_source`, `raw_count`, `library_size`,
`n_specimens`, `treatment`.

The recorded `timepoint` is a label, not an assumed numeric interval; include
an explicit order or elapsed time when needed. A patient-stage pseudobulk uses
summed raw counts divided by the summed full-gene library before applying the
declared transform. An average of cell-level log expression is a different
measure and must have a different `unit`. Neither is silently converted to the
other. Missing response remains missing; do not infer it from a sample name.
