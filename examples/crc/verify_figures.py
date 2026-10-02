#!/usr/bin/env python3
"""Check the entire frozen CRC figure inventory, PDFs and paired PNGs."""
import argparse
import csv
import json
from pathlib import Path

try:
    import pdfplumber
except ImportError as exc:
    raise SystemExit("PDF validation needs pdfplumber (pip install pdfplumber)") from exc


def main(root: Path) -> None:
    inventory = Path(__file__).with_name("figure_inventory.csv")
    with inventory.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 62:
        raise ValueError(f"figure inventory expected 62 PDF/PNG pairs, found {len(rows)}")
    report = []
    for row in rows:
        pdf = root / row["pdf"]
        png = root / row["png"]
        if not pdf.is_file() or not png.is_file():
            raise FileNotFoundError(f"missing pair: {pdf} / {png}")
        with pdfplumber.open(pdf) as document:
            if len(document.pages) != 1:
                raise ValueError(f"not one page: {pdf}")
            page = document.pages[0]
            if not page.chars:
                raise ValueError(f"no text glyphs: {pdf}")
            outside = sum(c["x0"] < -1 or c["x1"] > page.width + 1 or
                          c["top"] < -1 or c["bottom"] > page.height + 1
                          for c in page.chars)
            if outside:
                raise ValueError(f"{outside} glyphs outside page: {pdf}")
            report.append({"pdf": row["pdf"], "png": row["png"],
                           "generator": row["generator"],
                           "text_glyphs": len(page.chars), "outside_page": outside})
    output = root / "audit" / "figure_geometry_qa.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"VALID {len(report)} one-page PDF/PNG pairs; no outside-page text")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    main(parser.parse_args().work_root.resolve())
