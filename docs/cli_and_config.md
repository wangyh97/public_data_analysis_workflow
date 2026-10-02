# CLI and configuration

Paths that vary by computer or run belong on the command line (`--input`,
`--work-root`, `--source-dir`, `--out-dir`). Stable analysis decisions belong in
config: cell inclusion, normalization, gene sets, statistical tests and
multiple-testing families. Plot configs contain only appearance and layout.

The CRC runner resolves paths as explicit CLI value > existing work-root
`analysis_config.json` > executable discovery; its Python default is the
interpreter running the runner. Existing TCGA modules retain their own CLI and
environment-variable rules. No source-specific computer path should be committed. A run must save the
effective values, original command, code revision and input provenance under
its results directory. These values describe the run that actually happened,
not merely the defaults present in the repository.

Do not create a full config for each figure. One analysis config can produce
multiple related outputs; plot scripts read the resulting tables. A source
adapter should select `source` and study accessions through CLI or a short
study recipe. This is the extension rule for new adapters, not a claim that
the CRC case runner accepts arbitrary studies today. A special mapping file
is needed when deposited metadata cannot be interpreted using normal adapter
rules.

The existing TCGA modules retain their current CLI. The CRC case-study CLI is
documented in `examples/crc/README.md`; it stages frozen scripts into an
external work root and records each command in a run manifest.
