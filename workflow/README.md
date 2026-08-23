# workflow/

Minimal Snakemake proof-of-concept. **All logic lives in the R scripts**; the
Snakefile only wires dependencies, passes parameters, and defines outputs.

## Why so small?

The framework's CLI tools are designed as Snakemake-ready units:
deterministic inputs/outputs, no interactive steps, cwd-independent,
exit-code discipline. Wrapping them in rules is therefore trivial, and the
project never *depends* on Snakemake.

## Usage

```bash
export PDA_DATA_ROOT=/data/public_data
export PDA_RESULTS_ROOT=/project/results
snakemake -n --cores 1      # dry-run: verify the DAG resolves
snakemake --cores 1         # execute download -> prepare -> correlate -> plot
```

Change `configfile:` at the top of the Snakefile (or pass
`--configfile ...`) to run a different analysis. Cohorts to prepare are
enumerated from the analysis config's `datasets` list.

## Planned extensions (out of scope for now)

* per-rule `resources:` (mem/time) + a SLURM/cluster profile — resource
  requests belong here, never inside R scripts.
* content-addressed cache if metadata-based fingerprints ever prove insufficient.
