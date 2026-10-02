# CRC single-cell case study

This directory versions the code used for the existing CRC evidence set:
METTL3, STRAP and PTBP1 in GSE236581 and GSE205506, plus portal sensitivity
analyses and CellChat. The case code is fixed to the deposited study metadata
and is intentionally separate from source adapters intended for future GSEs.

## Inputs and entry point

The external `--work-root` holds `data/`, `derived/`, `audit/`, `tables/`,
`figures/` and `panel/`. The deposited CellResDB/scCT-DB files are supplied via
`--source-dir`; neither raw matrices nor local paths are committed. The runner
copies this versioned case code to the work root and generates
`analysis_config.json` using `analysis_parameters.json` plus CLI paths.

From the repository root in PowerShell:

```powershell
& 'C:/path/to/python.exe' examples/crc/run.py `
  --work-root 'C:/path/to/CRC_single_cell' `
  --source-dir 'C:/path/to/downloaded/single cell' `
  --rscript 'C:/path/to/Rscript.exe' `
  --node 'C:/path/to/node.exe' `
  --mode all
```

Modes: `stage` copies code only; `panels-plot` redraws the 27 current Fig1–3
PDFs/PNGs from panel source tables; `panels` regenerates tables from raw
matrices and draws them; `legacy-figures` regenerates the older summary and
sensitivity figures from prepared data; `legacy-full` starts with extraction;
`cellchat` prepares, infers and compares communication; `all` runs all those
phases; `validate` runs the panel numerical/PDF audit. Install the pinned
CellChat R runtime described in `docs/CellChat运行与设计.md` before running
`cellchat` or `all` on a new machine.

Add `--install-cellchat` to `cellchat` or `all` to install the case's pinned
runtime as part of the command. `--mode install-cellchat` performs just that
step. Add `--download-missing` to audit and fetch missing public input files
before full preparation. These steps require network access and are logged.

The `all` mode assumes the GEO/scCT archives and CellResDB files have already
been downloaded into the documented locations. Remote download/audit remains
available as frozen scripts `legacy/scripts/05_remote_audit.py` and
`06_download_missing.py`; run them explicitly if the raw input audit reports
missing files. The per-run `audit/runs/<timestamp>_<mode>/run_manifest.json`
contains the effective paths, parameters, Git revision, step commands, exit
codes and logs, including on failure.

To replay a run with its original committed code and arguments:

```powershell
python workflow/replay.py --manifest 'C:/path/to/audit/runs/<run>/run_manifest.json'
```

Use `--dry-run` to inspect it. Replay checks input identity and refuses
changed inputs unless `--allow-input-change` is supplied; that override marks
the rerun as non-exact. After plotting, the runner exports standardized v1
patient observation and statistical tables to `processed/contracts_v1/` and
validates their contracts.

## Figure groups and boundaries

- Legacy `figures/01`–`06`: `legacy/scripts/08_statistics_figures.R`.
- Legacy `figures/07`–`08`: `legacy/scripts/14_sensitivity.R`.
- Legacy `figures/09`–`11`: `legacy/scripts/22_cellchat_comparison.R`.
- Fig1 MHC-I and structural-gene scatter plots: panel 01 `script/plot.R`.
- Fig2 T-cell state heatmap: panel 02 `script/plot.R`.
- Fig3 trajectories, patient facets and paired boxplots: panel 03 `script/plot.R`.

The CellResDB and scCT-DB versions of a GEO study are alternate processing or
annotation versions of the same underlying patients, not independent cohorts.
GSE236581 response uses RECIST; GSE205506 uses pCR/non-pCR. These facts and
sample mappings are documented in `docs/` here, not in repository-wide docs.
See `figure_inventory.csv` for the exact expected figure filenames.
