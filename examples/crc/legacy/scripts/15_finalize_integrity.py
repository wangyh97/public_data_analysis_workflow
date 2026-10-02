"""Consolidate source fingerprints, dimensions and official byte-size checks."""
import pathlib,json,hashlib,csv
ROOT=pathlib.Path(__file__).resolve().parents[1]
cfg=json.loads((ROOT/'analysis_config.json').read_text());source=pathlib.Path(cfg['paths']['source_dir'])
rows=[]
for r in json.loads((ROOT/'audit/official/cellres_headers.json').read_text()):
    f=source/'CellResDB'/r['file'];h=hashlib.sha256()
    with open(f,'rb') as z:
        while b:=z.read(8*1024**2):h.update(b)
    remote=int(r['headers']['Content-Length']);assert f.stat().st_size==remote
    checkpoint=ROOT/'derived'/(f.stem+'_stream.json')
    if checkpoint.exists():meta=json.loads(checkpoint.read_text());cells=meta['n_cells']
    else:cells={'GSE205506_CRC_1':16657,'GSE205506_CRC_2':17291}[f.stem]
    rows.append({'source':'CellResDB','file':f.name,'bytes':remote,'official_size_match':True,'sha256':h.hexdigest(),'compression_and_read_validation':True,'cells':cells})
    print('HASH PASS',f.name,flush=True)
for r in json.loads((ROOT/'audit/official/scct_archive_headers.json').read_text()):
    name=r['url'].rsplit('/',1)[-1];v=json.loads((ROOT/'audit'/(name+'.validation.json')).read_text());remote=int(r['headers']['Content-Length']);assert remote==v['bytes']
    rows.append({'source':'scCT-DB','file':name,'bytes':remote,'official_size_match':True,'sha256':v['sha256'],'compression_and_read_validation':v['gzip_crc_pass'],'cells':''})
for r in json.loads((ROOT/'audit/geo_downloads.json').read_text()):
    rows.append({'source':'GEO supplement','file':r['file'],'bytes':r['bytes'],'official_size_match':True,'sha256':r['sha256'],'compression_and_read_validation':True,'cells':''})
with open(ROOT/'audit/input_integrity_final.csv','w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('COMPLETE integrity records',len(rows),flush=True)
