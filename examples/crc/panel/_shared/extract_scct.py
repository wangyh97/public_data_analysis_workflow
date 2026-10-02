"""Stream gene-by-cell CSV matrices without densifying the full dataset.
Check every row width, gene uniqueness, integer/nonnegative counts, barcode mapping.
Retain target/signature genes and full-library totals for reproducible normalization.
"""
import csv,gzip,pathlib,json,hashlib
import numpy as np
import sys
ROOT=pathlib.Path(sys.argv[1]).resolve()
GENES=json.loads((ROOT/'panel/_shared/extraction_genes.json').read_text())
out=ROOT/'derived/panel_v2';out.mkdir(exist_ok=True)
aud=[]
for folder in sorted((ROOT/'data/scCT-DB').iterdir()):
    if folder.name not in ['scCT_dataset_123-I','scCT_dataset_124-I','scCT_dataset_125-I']:continue
    mat=folder/'Matrix.csv.gz'
    if not mat.exists():continue
    prefix=out/folder.name
    if prefix.with_suffix('.json').exists():
        aud.append(json.loads(prefix.with_suffix('.json').read_text()))
        continue
    print('READ',folder.name,flush=True)
    with gzip.open(mat,'rt') as f:
        barcodes=next(csv.reader([f.readline()]))[1:];n=len(barcodes)
        lib=np.zeros(n);detected=np.zeros(n);mt=np.zeros(n);selected={};seen=set();integer=True;nonnegative=True
        for line in f:
            gene,vals=line.rstrip('\r\n').split(',',1);gene=gene.strip('"')
            if gene in seen:raise ValueError('Duplicate gene '+gene)
            seen.add(gene);v=np.fromstring(vals,sep=',')
            if len(v)!=n:raise ValueError('Row dimension mismatch '+gene)
            nonnegative &= bool(np.isfinite(v).all() and (v>=0).all())
            integer &= bool((np.abs(v-np.rint(v))<1e-8).all())
            lib+=v;detected+=v>0
            if gene.startswith('MT-'):mt+=v
            if gene in GENES:selected[gene]=v
    # CSVs are intentionally small and portable; one row per cell.
    with gzip.open(str(prefix)+'_expression.csv.gz','wt',newline='') as f:
        w=csv.writer(f);w.writerow(['barcode','library','nFeature','percent_mt']+list(selected))
        values=np.column_stack([lib,detected,100*mt/np.maximum(lib,1)]+list(selected.values()))
        w.writerows([b]+v.tolist() for b,v in zip(barcodes,values))
    info=list(csv.DictReader(open(folder/'Cell info.txt',encoding='utf-8-sig'),delimiter='\t'))
    ib=[r['Barcode'] for r in info]
    row={'dataset':folder.name,'n_cells':n,'n_genes':len(seen),'n_annotation_cells':len(ib),'barcode_set_equal':set(ib)==set(barcodes),'barcode_unique':len(set(barcodes))==n,'nonnegative_finite':nonnegative,'integer':integer,'missing_genes':sorted(set(GENES)-set(selected))}
    assert row['barcode_set_equal'] and row['barcode_unique'] and nonnegative and integer
    prefix.with_suffix('.json').write_text(json.dumps(row,indent=2))
    aud.append(row);print('DONE',row,flush=True)
(out/'scct_matrix_validation.json').write_text(json.dumps(aud,indent=2))
