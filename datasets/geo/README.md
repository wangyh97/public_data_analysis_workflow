# GEO source adapter

Source identity is GEO; a GSE accession, disease and treatment are selections
passed by CLI or study recipe. GEO contains many file formats, so parsing must
dispatch by deposited format and validate dimensions/metadata before emitting
the processing contract. The frozen GSE236581 reader used for the CRC figures
is retained under `examples/crc/panel/_shared/` as a study-specific reference.
This source-wide adapter is planned; the CRC reader is not silently applied to
arbitrary GSEs.
