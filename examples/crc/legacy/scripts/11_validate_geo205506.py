"""Validate all original GEO sample bundles, not just previously annotated cells.
Check gzip trailers for every member and matrix/barcode/feature dimensions.
The tar archive is retained; decompressed whole matrices are not written to disk.
"""
import tarfile,gzip,pathlib,csv,json,collections
ROOT=pathlib.Path(__file__).resolve().parents[1];rows=[];samples={}
with tarfile.open(ROOT/'data/GEO/GSE205506_RAW.tar') as t:
    for m in t:
        if not m.isfile():continue
        sample=m.name.split('_')[0];rec=samples.setdefault(sample,{'sample':sample})
        with gzip.GzipFile(fileobj=t.extractfile(m)) as f:
            if '_matrix.mtx' in m.name:
                line=f.readline()
                assert line.startswith(b'%%MatrixMarket matrix coordinate')
                line=f.readline()
                while line.startswith(b'%'):line=f.readline()
                ng,nc,nnz=map(int,line.split());rec.update(matrix_genes=ng,matrix_cells=nc,nnz=nnz)
                lines=0
                while b:=f.read(8*1024**2):lines+=b.count(b'\n')
                if lines!=nnz:raise ValueError((sample,lines,nnz))
            else:
                n=0
                while b:=f.read(8*1024**2):n+=b.count(b'\n')
                rec['features' if '_features.tsv' in m.name else 'barcodes']=n
        rows.append({'member':m.name,'bytes':m.size,'gzip_crc_pass':True})
        print('PASS',m.name,flush=True)
for s in samples.values():
    s['dimensions_match']=s['features']==s['matrix_genes'] and s['barcodes']==s['matrix_cells']
    assert s['dimensions_match']
with open(ROOT/'audit/GSE205506_full_sample_validation.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(next(iter(samples.values()))));w.writeheader();w.writerows(samples.values())
(ROOT/'audit/GSE205506_full_member_validation.json').write_text(json.dumps(rows,indent=2))
print('COMPLETE samples',len(samples),'cells',sum(s['barcodes'] for s in samples.values()),flush=True)
