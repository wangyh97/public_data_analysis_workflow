# Add a plot

Create `plotting/<type>/` and document required input table fields. The CLI
accepts input table(s), output directory and an optional stable style config.
It reads stored estimates and test results without running tests. Keep
PDF/PNG export and figure size together in the plotting module; write the
effective style and exact inputs to the run manifest. A single plot type can
produce multiple genes, cohorts or formats from one invocation.
