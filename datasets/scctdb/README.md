# scCT-DB source adapter

scCT-DB is a data provider. Portal dataset IDs and underlying GEO study IDs
are separate fields. The frozen CRC matrix reader and sample mapping are
versioned under `examples/crc/panel/_shared/`. A future source-wide adapter
must accept selected portal IDs at runtime, validate barcodes and raw integer
counts, and emit the documented processing contracts. Portal reannotation of
a GEO study is not an independent biological cohort.
