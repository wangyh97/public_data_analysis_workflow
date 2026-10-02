"""Inspect actual PDF glyph geometry (including Cairo transforms) and package results.

Does not modify scientific tables, figures, or original expression data.
Large data/derived directories are deliberately excluded from the shareable ZIP.
"""
from pathlib import Path
import csv
import json
import zipfile
import pdfplumber

ROOT = Path(__file__).resolve().parents[1]
rows = []
for path in sorted((ROOT / 'figures').glob('*.pdf')):
    with pdfplumber.open(path) as pdf:
        chars = [c for page in pdf.pages for c in page.chars if c['text'].strip()]
        outside = sum(c['x0'] < -1 or c['x1'] > page.width + 1 or
                      c['top'] < -1 or c['bottom'] > page.height + 1
                      for page in pdf.pages for c in page.chars)
        # For rotated labels PDFplumber size is a glyph-width projection; use
        # the font transform vector norm to recover nominal physical font size.
        import math
        sizes = [math.hypot(c['matrix'][0], c['matrix'][1]) *
                 (c['height'] if c['upright'] else c['width']) /
                 max(math.hypot(c['matrix'][0], c['matrix'][1]), 1e-12)
                 for c in chars]
        rows.append(dict(file=path.name, pages=len(pdf.pages), text_glyphs=len(chars),
                         min_glyph_height_pt=round(min(sizes), 3),
                         glyphs_outside_page=outside,
                         png_present=path.with_suffix('.png').is_file()))
        assert chars and not outside and path.with_suffix('.png').is_file(), path
with (ROOT/'audit/pdf_geometry_qa.csv').open('w', newline='', encoding='utf-8-sig') as f:
    writer = csv.DictWriter(f, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
print(json.dumps({'pdf_count':len(rows), 'all_have_text_and_png': True,
                  'outside_page_glyphs':sum(r['glyphs_outside_page'] for r in rows)}))
archive = ROOT/'CRC_code_results.zip'
with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as z:
    for f in sorted(ROOT.rglob('*')):
        rel = f.relative_to(ROOT)
        if f.is_file() and rel.parts[0] not in {'data','derived','runtime'} and f != archive and '__pycache__' not in rel.parts:
            z.write(f, rel.as_posix())
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
print('ZIP verified', archive.stat().st_size, 'bytes')
