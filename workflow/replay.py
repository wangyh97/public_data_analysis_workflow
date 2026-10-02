#!/usr/bin/env python3
"""Replay a CRC run manifest with the exact committed case-study code.

Usage: python workflow/replay.py --manifest /path/to/run_manifest.json
The original working data directory is reused. Add --dry-run to inspect the
resolved historical code and command without changing files.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--allow-input-change", action="store_true",
                   help="run despite changed input identity; marked non-exact")
    args = p.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    commit = manifest.get("git_commit")
    if not commit or not all(c in "0123456789abcdef" for c in commit) or len(commit) != 40:
        p.error("manifest lacks a full Git commit; exact code replay is unavailable")
    work_root = Path(manifest["paths"]["work_root"])
    changes = []
    for rel, expected in manifest.get("input_identity", {}).items():
        path = work_root / rel
        if not path.is_file() or digest(path) != expected["sha256"]:
            changes.append(rel)
    if changes and not args.allow_input_change:
        p.error("input identity changed: " + ", ".join(changes))
    command = manifest["command"]
    original_script = Path(command[0])
    if original_script.name != "run.py" or "examples/crc" not in original_script.as_posix():
        p.error("manifest is not from examples/crc/run.py")
    git = ["git", "-c", f"safe.directory={REPO.as_posix()}", "-C", str(REPO)]
    subprocess.run(git + ["cat-file", "-e", f"{commit}^{{commit}}"], check=True)
    replay_args = command[1:]
    # Old run names must never overwrite a historical manifest.
    if "--run-name" in replay_args:
        ix = replay_args.index("--run-name")
        del replay_args[ix:ix + 2]
    # Resolve all paths from the historical effective manifest, not from the
    # possibly changed config currently present in the work directory.
    saved = manifest["paths"]
    for flag, key in (("--work-root", "work_root"), ("--source-dir", "source_dir"),
                      ("--rscript", "rscript"), ("--python", "python"),
                      ("--node", "node")):
        value = saved.get(key)
        if value:
            if flag in replay_args:
                ix = replay_args.index(flag)
                replay_args[ix + 1] = value
            else:
                replay_args.extend((flag, value))
    print("Git revision:", commit)
    print("Work root:", work_root)
    print("Input changes:", changes or "none")
    print("Replay arguments:", replay_args)
    if args.dry_run:
        return
    archive = subprocess.check_output(git + ["archive", "--format=tar", commit,
                                           "examples/crc", "processing/contracts"])
    with tempfile.TemporaryDirectory(prefix="pda_crc_replay_") as temporary:
        checkout = Path(temporary)
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as tar:
            for member in tar:
                target = (checkout / member.name).resolve()
                if not target.is_relative_to(checkout) or member.issym() or member.islnk():
                    raise RuntimeError("unsafe archived path: " + member.name)
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with tar.extractfile(member) as src, target.open("wb") as dst:
                        dst.write(src.read())
        archived_runner = checkout / "examples" / "crc" / "run.py"
        environment = os.environ.copy()
        environment["PDA_REPLAY_GIT_COMMIT"] = commit
        environment["PDA_REPLAY_SOURCE_MANIFEST"] = str(args.manifest.resolve())
        if changes:
            environment["PDA_REPLAY_INPUT_CHANGED"] = "true"
        subprocess.run([sys.executable, str(archived_runner), *replay_args],
                       env=environment, check=True)


if __name__ == "__main__":
    main()
