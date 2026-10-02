#!/usr/bin/env python3
"""Validate version-1 long-table contracts without optional dependencies."""
import argparse
import csv
import math
from pathlib import Path

CONTRACTS = {
    "observations": {
        "required": ("source", "study_id", "patient_id", "timepoint", "compartment",
                     "feature", "value", "unit", "n_cells"),
        "key": ("source", "study_id", "patient_id", "timepoint", "compartment", "feature"),
    },
    "statistics": {
        "required": ("analysis_id", "source", "study_id", "comparison", "feature",
                     "method", "n", "estimate", "p_value", "adjusted_p",
                     "adjustment_family", "status"),
        "key": ("analysis_id", "source", "study_id", "comparison", "feature",
                "adjustment_family"),
    },
}


def validate(path: Path, kind: str) -> int:
    schema = CONTRACTS[kind]
    with path.open(newline="", encoding="utf-8-sig") as stream:
        sample = stream.read(4096)
        stream.seek(0)
        dialect = csv.excel_tab if path.suffix.lower() == ".tsv" else csv.excel
        reader = csv.DictReader(stream, dialect=dialect)
        if reader.fieldnames is None:
            raise ValueError("empty table")
        missing = set(schema["required"]) - set(reader.fieldnames)
        if missing:
            raise ValueError("missing columns: " + ", ".join(sorted(missing)))
        seen = set()
        count = 0
        for line, row in enumerate(reader, 2):
            if any(not row[col] for col in schema["key"]):
                raise ValueError(f"line {line}: empty key")
            key = tuple(row[col] for col in schema["key"])
            if key in seen:
                raise ValueError(f"line {line}: duplicate key {key}")
            seen.add(key)
            nfield = "n_cells" if kind == "observations" else "n"
            n = int(row[nfield])
            if n < 0:
                raise ValueError(f"line {line}: negative {nfield}")
            if kind == "observations" and not math.isfinite(float(row["value"])):
                raise ValueError(f"line {line}: non-finite value")
            if kind == "statistics" and row["status"] == "tested":
                for name in ("estimate", "p_value", "adjusted_p"):
                    v = float(row[name])
                    if not math.isfinite(v):
                        raise ValueError(f"line {line}: non-finite {name}")
                    if name != "estimate" and not 0 <= v <= 1:
                        raise ValueError(f"line {line}: {name} outside [0,1]")
            count += 1
    if not count:
        raise ValueError("table has no data rows")
    return count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", choices=sorted(CONTRACTS), required=True)
    parser.add_argument("--input", type=Path, required=True)
    args = parser.parse_args()
    print(f"VALID {args.contract}: {validate(args.input, args.contract)} rows")
