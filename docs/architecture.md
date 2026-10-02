# Architecture

The intended module boundary is:

```
source adapter -> processing/contract -> analysis table -> plotting -> figure
```

- `datasets/<source>/` obtains and parses a source format. A source is GEO,
  CellResDB, scCT-DB, TCGA, etc.; CRC or a GSE number is a study selection.
- `processing/` harmonizes identifiers, annotations and expression units. It
  emits one of the versioned contracts in `data_contracts/`.
- `analyses/` consumes validated processed data and writes numerical tables.
  The patient or specimen, never an individual cell, is the independent unit
  when the analysis is patient-level.
- `plotting/` consumes prepared data and statistics. Changing only style does
  not rerun data processing or hypothesis tests.

The existing TCGA adapter uses these general modules. The CRC example keeps
its validated historical scripts under `examples/crc/legacy/`, and its three
new panels separate data preparation from plotting. It exports two v1 contract
tables so future analyses can consume the same patient-level interface.
Source-wide GEO, CellResDB and scCT-DB adapters are still planned; selecting a
different accession in the CRC case runner does not yet generalize it.

`workflow/` holds orchestration and historical replay; `utils/` holds existing
common helpers. Large raw data, caches and figures live outside Git.
