# Add an analysis

Create `analyses/<method>/` with a noninteractive CLI and a short README.
Read a named data contract, validate its columns/units, state the independent
unit and inclusion rules, and write one row per hypothesis plus all plot-ready
values to tables. Record raw p values and the exact adjustment family.
Analysis code must not read plotting colors or draw figures. Add a numerical
regression check with a small known dataset and include the module in the
workflow dependency graph.
