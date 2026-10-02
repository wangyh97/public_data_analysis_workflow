"""Stream the complete GEO Matrix Market matrix into selected-gene counts.
All entries contribute to library totals; all raw cell barcodes are checked.
No sampling of matrix entries or cells. Peak arrays are bounded and disk-backed.
"""
import pathlib,gzip,csv,json,re,time
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1];src=ROOT/'data/GEO'
GENES=re.findall(r"'([^']+)'",(ROOT/'scripts/02_extract_cellres.R').read_text().split('genes <- c(')[1].split(')\n')[0])
prefix=ROOT/'derived/GEO_GSE236581';scratch=ROOT/'data/stream_scratch/GEO_GSE236581.bin'
features=[line.rstrip('\n').split('\t') for line in gzip.open(src/'GSE236581_features.tsv.gz','rt')]
symbols=[r[1] if len(r)>1 else r[0] for r in features]
bars=[line.strip().strip('"') for line in gzip.open(src/'GSE236581_barcodes.tsv.gz','rt')]
present=[g for g in GENES if g in symbols];lookup=np.full(len(symbols),-1,dtype=np.int32)
for j,g in enumerate(present):
    ix=[i for i,s in enumerate(symbols) if s==g]
    if len(ix)!=1:raise ValueError('Ambiguous duplicate gene '+g)
    lookup[ix[0]]=j
is_mt=np.array([g.startswith('MT-') for g in symbols]);n=len(bars)
e=np.memmap(scratch,dtype='float64',mode='w+',shape=(n,len(present)));e[:]=0
lib=np.zeros(n);det=np.zeros(n,dtype='int32');mt=np.zeros(n);seen=0
with gzip.open(src/'GSE236581_counts.mtx.gz','rb') as f:
    line=f.readline()
    if not line.startswith(b'%%MatrixMarket matrix coordinate'):raise ValueError('Unexpected matrix format')
    line=f.readline()
    while line.startswith(b'%'):line=f.readline()
    ng,nc,nnz=map(int,line.split());assert ng==len(symbols) and nc==n
    print('Dimensions',ng,nc,nnz,flush=True)
    carry=b'';last_col=-1;last_row=-1
    while block:=f.read(16*1024**2):
        block=carry+block;end=block.rfind(b'\n');carry=block[end+1:]
        v=np.fromstring(block[:end].decode('ascii'),sep=' ',dtype=np.int64).reshape(-1,3)
        r=v[:,0]-1;c=v[:,1]-1;x=v[:,2]
        assert (r>=0).all() and (r<ng).all() and (c>=0).all() and (c<nc).all() and np.isfinite(x).all() and (x>=0).all() and (x==np.rint(x)).all()
        # The deposited matrix is CSC-ordered. Enforce ordering/coordinate uniqueness
        # before using fast segmented reductions instead of slow scattered addition.
        assert (np.diff(c)>=0).all() and ((np.diff(c)>0)|(np.diff(r)>0)).all()
        assert c[0]>last_col or (c[0]==last_col and r[0]>last_row)
        last_col=int(c[-1]);last_row=int(r[-1])
        starts=np.r_[0,np.flatnonzero(np.diff(c))+1];cols=c[starts]
        lib[cols]+=np.add.reduceat(x,starts);det[cols]+=np.add.reduceat((x>0).astype('int32'),starts)
        mt[cols]+=np.add.reduceat(x*is_mt[r],starts)
        j=lookup[r];mask=j>=0;e[c[mask],j[mask]]+=x[mask];seen+=len(v)
        if seen%50000000<len(v):e.flush();print('Entries',seen,'/',nnz,flush=True)
    assert not carry.strip() and seen==nnz
e.flush()
with gzip.open(str(prefix)+'_expression.csv.gz','wt',newline='') as f:
    w=csv.writer(f);w.writerow(['barcode','library','nFeature','percent_mt']+present)
    for j,b in enumerate(bars):w.writerow([b,lib[j],det[j],100*mt[j]/max(lib[j],1)]+e[j].tolist())
summary={'n_genes':ng,'n_cells':nc,'nnz':nnz,'entries_read':seen,'gzip_crc_pass':True,'genes_missing':sorted(set(GENES)-set(present))}
(ROOT/'audit/GEO_GSE236581_matrix_validation.json').write_text(json.dumps(summary,indent=2))
e._mmap.close();del e;scratch.unlink();print('DONE',summary,flush=True)
