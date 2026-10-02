# GEO_SC_CRC adapter

This adapter connects the CRC single-cell panel pipeline to the workflow
without duplicating the large GEO matrices. `panel_root` points to the
existing `outputs/CRC_single_cell/panel` directory. Cell filtering,
pseudobulk aggregation, response definitions, and statistics remain in the
validated panel scripts. The wrapper stages source tables, figures, logs, and
provenance below `PDA_RESULTS_ROOT`.

```powershell
$env:PDA_RESULTS_ROOT='results'
$env:CRC_PANEL_ROOT='C:/path/to/outputs/CRC_single_cell/panel'
Rscript analyses/crc_single_cell/run.R --config config/analyses/crc_single_cell.yaml
```
