"""Extract all genes for selected scCT cells to the same sparse triplet format."""
from pathlib import Path
import csv,gzip,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];out=ROOT/'data/CellChat'
sel=pd.read_csv(out/'selected_cells.csv')
for ds,m in sel[sel.dataset!='GEO'].groupby('dataset'):
 states={};genes=[]
 with gzip.open(ROOT/'data/scCT-DB'/ds/'Matrix.csv.gz','rt') as f:
  bars=next(csv.reader([f.readline()]))[1:];lookup={b:i for i,b in enumerate(bars)}
  for job,g in m.groupby('job'):
   states[job]=dict(file=open(out/(job+'_counts.bin'),'wb'),ix=np.array([lookup[b] for b in g.barcode]),cols=g.local_col.to_numpy(),nnz=0,n=len(g))
  for line in f:
   gene,values=line.rstrip('\r\n').split(',',1);genes.append(gene.strip('"'));v=np.fromstring(values,sep=',')
   assert len(v)==len(bars) and np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all()
   for s in states.values():
    a=v[s['ix']];nz=a>0
    trip=np.column_stack([np.repeat(len(genes),nz.sum()),s['cols'][nz],a[nz]]).astype('<u4')
    s['file'].write(trip.tobytes());s['nnz']+=int(nz.sum())
 for job,s in states.items():
  s['file'].close();(out/(job+'_genes.txt')).write_text('\n'.join(genes)+'\n')
  (out/(job+'_extraction.json')).write_text(json.dumps(dict(job=job,genes=len(genes),nnz=s['nnz'],cells=s['n'])))
 print('COMPLETE',ds,len(states),'samples',flush=True)
