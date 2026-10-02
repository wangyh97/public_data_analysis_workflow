# CellResDB source adapter

CellResDB is a data provider. Study accession and tumor type are selections,
not part of the adapter name. The frozen CRC Seurat-RDS extraction and
annotation audit are versioned under `examples/crc/legacy/scripts/` for exact
reproduction. A future source-wide adapter must expose its own dataset-ID
selection and validate the deposited RDS structure before applying common
processing contracts.
