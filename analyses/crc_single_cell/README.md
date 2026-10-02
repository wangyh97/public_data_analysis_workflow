# CRC single-cell analysis module

This integration layer follows the workflow CLI, output, and provenance
conventions while reusing the validated panel scripts. Results are staged at
`$PDA_RESULTS_ROOT/crc_single_cell_panels/panels/` with one directory for each
panel and its `source/`, `script/`, `figure/`, and README files. Large raw GEO
matrices remain outside the repository.
