# tests/ (deferred)

Per project decision, automated tests are **deferred** for this phase.
Planned smoke test (to run on HPC once enabled):

1. `prepare.R --cohorts SKCM` on a small cohort; verify processed contract files.
2. Missing-gene behaviour: add FAKEGENE to an analysis config -> warn, not fail.
3. Config-change rerun matrix:
   - edit plotting YAML only  -> only plot re-runs
   - edit analysis YAML       -> analyze + plot re-run
   - unchanged everything     -> download/prepare reuse via fingerprints
4. `snakemake -n` DAG matches expectation.

A synthetic mini expression/metadata fixture will live here so steps 2-4 can
run fully offline without the 1.3 GB bulk file.
