# Run a new study

1. Identify source (`geo`, `cellresdb`, `scctdb`, `tcga`, etc.) and the accession
   or dataset IDs. Keep these two concepts separate; one study can appear in
   multiple portals.
2. Provide raw/source and output paths via CLI. Use the source adapter to
   download/parse and the processing module to generate a validated contract.
3. Choose an analysis module and a stable rule config. Specify the biological
   independent unit, grouping, response endpoint and multiple-testing family.
4. Plot the resulting tables with a plotting module. Change style without
   recalculating the analysis.
5. Inspect the run manifest and audit tables before interpreting figures.

The frozen CRC study at `examples/crc/` demonstrates historical figure
reproduction. It is not a source adapter: changing a GSE number in that case
does not magically reinterpret the new study's metadata. New studies should
use the generic modules and add a small mapping only for genuine exceptions.
