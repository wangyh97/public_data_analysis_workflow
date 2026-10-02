# Expression matrix contract (v1)

An expression matrix has features as rows and independent specimens as
columns. Column identifiers match `sample_metadata.sample_id` exactly and are
unique. Row identifiers match `feature_metadata.feature_id` and are unique
after the declared duplicate-resolution rule. The manifest states whether
values are raw counts, TPM, log TPM, or another unit; scripts never infer a
unit from numeric range alone.

For single-cell experiments, the raw gene-by-cell matrix remains a source
artifact. Patient-level analyses first generate the observation contract,
recording the cell compartment, cell count and aggregation rule. Distinct
normalizations (for example mean log1p CP10K versus log2 pseudobulk CPM+1)
must not share an ambiguous `expression` unit label.
