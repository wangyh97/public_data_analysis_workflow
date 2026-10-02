"""Memory-bounded reader for the XDR v3 Seurat RDS serialization used here.
The complete stream is parsed and gzip CRC checked, but large redundant assays
are discarded. Sparse RNA counts are reduced using on-disk row-index/target arrays.
This is deliberately format-specific; unknown R types fail loudly, never guessed.
Validation against native R readRDS on GSE205506 is required before using large files.
"""
import pathlib,gzip,lzma,struct,json,csv,sys,os,hashlib
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
SOURCE=pathlib.Path(json.loads((ROOT/'analysis_config.json').read_text())['paths']['source_dir'])/'CellResDB'
import re
GENES=re.findall(r"'([^']+)'",(ROOT/'scripts/02_extract_cellres.R').read_text().split('genes <- c(')[1].split(')\n')[0])
SCRATCH=ROOT/'data'/'stream_scratch';SCRATCH.mkdir(exist_ok=True)
class MetadataComplete(Exception):pass
class Parser:
    def __init__(self,f,label):self.f=f;self.refs=[];self.label=label;self.counts={};self.metadata=None;self.selected=None;self.integer=True
    def read(self,n):
        b=self.f.read(n)
        if len(b)!=n:raise EOFError((n,len(b),self.f.tell()))
        return b
    def i(self):return struct.unpack('>i',self.read(4))[0]
    def skip(self,n):
        while n: k=min(n,8*1024**2);self.read(k);n-=k
    def length(self):
        n=self.i()
        return ((self.i()<<32)|self.i()) if n==-1 else n
    def numeric(self,t,n,path):
        dtype='>i4' if t in (10,13) else '>f8';size=4 if t in (10,13) else 8
        iscounts='/counts/' in path and '/assays/' in path
        name=path.rsplit('/',1)[-1]
        if iscounts and name=='i':
            p=SCRATCH/(self.label+'.indices.bin');a=np.memmap(p,dtype='int32',mode='w+',shape=n)
            for s in range(0,n,1000000):k=min(1000000,n-s);a[s:s+k]=np.frombuffer(self.read(k*size),dtype=dtype)
            a.flush();self.indices=a;self.indexpath=p;return {'disk_indices':n}
        if iscounts and name=='x':
            c=self.counts;genes=c['Dimnames'][0];bars=c['Dimnames'][1];p=np.asarray(c['p'],dtype=np.int64);nc=len(bars)
            assert n==len(self.indices) and p[-1]==n and len(p)==nc+1
            present=[g for g in GENES if g in genes];lookup=np.full(len(genes),-1,dtype=np.int32)
            for j,g in enumerate(present):lookup[genes.index(g)]=j
            ep=SCRATCH/(self.label+'.expression.bin');e=np.memmap(ep,dtype='float64',mode='w+',shape=(nc,len(present)));e[:]=0
            lib=np.zeros(nc);det=np.diff(p);mt=np.zeros(nc);is_mt=np.array([g.startswith('MT-') for g in genes])
            for s in range(0,n,500000):
                k=min(500000,n-s);v=np.frombuffer(self.read(k*8),dtype='>f8').astype('float64');ix=self.indices[s:s+k]
                if not np.isfinite(v).all() or (v<0).any():raise ValueError('Invalid counts')
                self.integer &= bool((np.abs(v-np.rint(v))<1e-8).all())
                cols=np.searchsorted(p,np.arange(s,s+k),side='right')-1
                np.add.at(lib,cols,v)
                mask=is_mt[ix];np.add.at(mt,cols[mask],v[mask])
                j=lookup[ix];mask=j>=0;e[cols[mask],j[mask]]=v[mask]
            e.flush();del ix;self.indices._mmap.close();del self.indices;self.indexpath.unlink()
            self.selected=(present,bars,e,lib,det,mt,ep)
            print('COUNTS',self.label,nc,len(genes),n,flush=True)
            return {'reduced':True,'nnz':n}
        if n>2000000 or ('/data/' in path or '/scale.data' in path or '/graphs/' in path):self.skip(n*size);return {'skipped':n}
        a=np.frombuffer(self.read(n*size),dtype=dtype).astype('int64' if size==4 else 'float64')
        return a
    def obj(self,path='root'):
        flags=self.i();t=flags&255;attr=bool(flags&512);tag=bool(flags&1024)
        if t==254:return None
        if t==255:
            idx=flags>>8
            if not idx:idx=self.i()
            return self.refs[idx-1]
        if t==1:
            v=self.obj(path+'/symbol');self.refs.append(v);return v
        if t==9:
            n=self.i();return None if n==-1 else self.read(n).decode('utf-8',errors='replace')
        if t in (2,6):
            at=self.obj(path+'/attrs') if attr else None
            key=self.obj(path+'/tag') if tag else None
            pp=path+'/'+str(key) if key else path
            v=self.obj(pp)
            if '/assays/' in pp and '/counts/' in pp:self.counts[str(key)]=v
            if key=='meta.data':
                self.metadata=v
                raise MetadataComplete()
            tail=self.obj(path)
            pairs=[(key,v)]+(tail if isinstance(tail,list) else [])
            return pairs
        if t in (10,13,14):v=self.numeric(t,self.length(),path)
        elif t in (16,19,20):v=[self.obj(path+'/'+str(j)) for j in range(self.length())]
        elif t==24:
            n=self.length();self.skip(n);v={'raw_bytes':n}
        elif t==25:v={'S4':True}
        elif t in (241,242):
            v={'namespace':self.obj(path+'/namespace')};self.refs.append(v);return v
        elif t in (253,252,251):return {'special':t}
        else:raise ValueError(('Unsupported R type',t,flags,path,self.f.tell()))
        if attr:
            at=self.obj(path+'/attributes')
            return {'value':v,'attrs':dict(at or [])}
        return v
def values(x):
    if isinstance(x,dict) and 'value' in x:
        a=x['attrs'];v=x['value']
        if 'levels' in a:
            lev=a['levels'];return [lev[int(i)-1] if i>0 else '' for i in v]
        return values(v)
    if isinstance(x,np.ndarray):return x.tolist()
    return x
def run(path):
    label=path.stem;prefix=ROOT/'derived'/label
    if pathlib.Path(str(prefix)+'_stream.json').exists():return
    print('READ',label,flush=True)
    opener=gzip.open if path.open('rb').read(2)==b'\x1f\x8b' else lzma.open
    with opener(path,'rb') as f:
        assert f.read(2)==b'X\n';p=Parser(f,label);version=p.i();p.i();p.i();assert version in (2,3)
        if version==3:p.read(p.i())
        try:p.obj()
        except MetadataComplete:pass
        while f.read(8*1024**2):pass
    m=p.metadata;assert m is not None
    attrs=m['attrs'];names=attrs['names'];cols=[values(v) for v in m['value']];bars=values(attrs['row.names'])
    g,cb,e,lib,det,mt,ep=p.selected
    assert len(bars)==len(cb) and bars==cb
    with gzip.open(str(prefix)+'_metadata.csv.gz','wt',newline='') as f:
        w=csv.writer(f);w.writerow(['barcode']+names);w.writerows([b]+list(row) for b,row in zip(bars,zip(*cols)))
    with gzip.open(str(prefix)+'_expression.csv.gz','wt',newline='') as f:
        w=csv.writer(f);w.writerow(['barcode','library','nFeature','percent_mt']+g)
        for j,b in enumerate(cb):w.writerow([b,lib[j],det[j],100*mt[j]/max(lib[j],1)]+e[j].tolist())
    summary={'dataset':label,'n_cells':len(cb),'n_genes':len(p.counts['Dimnames'][0]),'integer':p.integer,'gzip_crc_pass':True,'barcode_order_equal':True,'genes_missing':sorted(set(GENES)-set(g))}
    pathlib.Path(str(prefix)+'_stream.json').write_text(json.dumps(summary,indent=2))
    e._mmap.close();del e;p.selected=None;ep.unlink();print('DONE',summary,flush=True)
if __name__=='__main__':
    files=sorted(SOURCE.glob(sys.argv[1] if len(sys.argv)>1 else '*Colorectal*.rds'))
    for path in files:run(path)
