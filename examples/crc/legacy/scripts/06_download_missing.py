"""Download missing original processed data and annotations, atomically and resumably.
Does not overwrite original user files. 1 MiB chunks; validates Content-Length and
records SHA256. Network access is required. No controlled raw FASTQ is requested.
"""
import pathlib,urllib.request,hashlib,json,time,gzip,tarfile
ROOT=pathlib.Path(__file__).resolve().parents[1];out=ROOT/'data/GEO';out.mkdir(exist_ok=True)
urls=['https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205506/suppl/GSE205506_RAW.tar',
      'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205506/suppl/filelist.txt']
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236581/suppl/'
urls += [base+'GSE236581_'+f for f in ['CRC-ICB_metadata.txt.gz','barcodes.tsv.gz','features.tsv.gz','counts.mtx.gz','VDJ_merge.txt.gz']]
records=[]
for url in urls:
    name=url.rsplit('/',1)[-1];dest=out/name;part=out/(name+'.part')
    with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=90) as r:expected=int(r.headers['Content-Length'])
    for attempt in range(8):
        if dest.exists() and dest.stat().st_size==expected:break
        start=part.stat().st_size if part.exists() else 0
        try:
            req=urllib.request.Request(url,headers={'Range':f'bytes={start}-'} if start else {})
            with urllib.request.urlopen(req,timeout=120) as r:
                if start and r.status!=206:raise ValueError('Server did not honor resume Range')
                with open(part,'ab' if start else 'wb') as f:
                    while b:=r.read(1024**2):f.write(b)
            assert part.stat().st_size==expected
            part.replace(dest);break
        except Exception as e:print(name,'retry',attempt,str(e),flush=True);time.sleep(3)
    if not dest.exists():raise RuntimeError('Download incomplete '+name)
    h=hashlib.sha256()
    with open(dest,'rb') as f:
        while b:=f.read(8*1024**2):h.update(b)
    if name.endswith('.gz'):
        with gzip.open(dest,'rb') as f:
            while f.read(8*1024**2):pass
    if name.endswith('.tar'):
        with tarfile.open(dest) as t:
            members=[{'name':m.name,'bytes':m.size} for m in t if m.isfile()]
        (ROOT/'audit/GSE205506_raw_members.json').write_text(json.dumps(members,indent=2))
    row={'file':name,'url':url,'bytes':expected,'sha256':h.hexdigest(),'status':'downloaded_or_verified'};records.append(row)
    (ROOT/'audit/geo_downloads.json').write_text(json.dumps(records,indent=2));print('DONE',row,flush=True)
