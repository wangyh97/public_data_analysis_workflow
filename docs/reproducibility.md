# Reproducibility

The CRC runner writes a manifest beside outputs with absolute work/source
paths, the exact command, effective parameters, Git revision, SHA-256 hashes
of selected small source manifests and audit tables, executed steps and log
locations. A failure after manifest creation records its status and the
failing step's log. Large raw matrices remain outside Git; the current runner
does not hash every matrix during each run.

Verification has three levels: input integrity (file IDs/cell/sample counts),
numerical results (patient counts, estimates, p and q values), and visual
outputs (expected PDF/PNG set, readable text, no clipping). A workflow is
called reproduced only after all three pass for the selected scope. The CRC
case inventory is in `examples/crc/README.md`; its original data and results
remain in the external work root.

For a case run recorded by `examples/crc/run.py`, replay the historical
instruction with one command:

```text
python workflow/replay.py --manifest /path/to/audit/runs/<run>/run_manifest.json
```

Replay checks recorded small input manifests/metadata by SHA-256, extracts
the exact recorded Git commit into a temporary directory, and executes the
stored arguments. It writes a new run manifest linked to the original.
`--dry-run` shows what would execute. Changed inputs are refused by default;
`--allow-input-change` is marked as a non-exact rerun. Historical case scripts
write to the same work root, so archive previous figures/tables you need to
retain before replaying.
