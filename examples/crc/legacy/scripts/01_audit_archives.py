"""Audit local scCT-DB archives, safely extract CRC members, save source manifests.
Uses only Python standard library; original downloads are read-only.
Run from the analysis output directory: python scripts/01_audit_archives.py
"""
import csv, io, json, tarfile, gzip, hashlib, pathlib, shutil, urllib.request, concurrent.futures, datetime
ROOT = pathlib.Path(__file__).resolve().parents[1]
SOURCE = pathlib.Path(json.loads((ROOT/'analysis_config.json').read_text())['paths']['source_dir'])
AUDIT=ROOT/'audit'; DATA=ROOT/'data'/'scCT-DB'
AUDIT.mkdir(exist_ok=True); DATA.mkdir(parents=True,exist_ok=True)
def save_rows(path, rows):
    if rows:
        with open(path,'w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
rows=[]; selected=set(); inventory=[]
with tarfile.open(SOURCE/'scCT-DB'/'SampleInfo.tar.gz') as t:
    for m in t:
        if not m.isfile(): continue
        rr=list(csv.DictReader(io.StringIO(t.extractfile(m).read().decode('utf-8-sig')),delimiter='\t'))
        crc=[r for r in rr if any(k in r.get('Cancer subtype','').lower() for k in ['colorectal','colon cancer','rectal'])]
        if crc:
            selected.add(m.name.split('/')[0]); rows.extend(crc)
save_rows(AUDIT/'scct_crc_samples.csv',rows)
print('CRC dataset IDs:',sorted(selected),flush=True)
for archive in ['SampleInfo.tar.gz','CellInfo.tar.gz','Matrix.tar.gz']:
    path=SOURCE/'scCT-DB'/archive
    print('Scanning',archive,flush=True)
    # Read gzip to EOF even after tar end markers; this verifies gzip CRC/trailer.
    h=hashlib.sha256(); total=0
    with open(path,'rb') as f:
        while b:=f.read(8*1024*1024): h.update(b)
    with gzip.open(path,'rb') as gz:
        with tarfile.open(fileobj=gz,mode='r|') as t:
            for m in t:
                inventory.append({'archive':archive,'member':m.name,'bytes':m.size,'crc_selected':m.name.split('/')[0] in selected})
                if m.isfile() and m.name.split('/')[0] in selected:
                    dest=(DATA/m.name).resolve()
                    if not dest.is_relative_to(DATA.resolve()): raise ValueError('Unsafe tar member')
                    dest.parent.mkdir(parents=True,exist_ok=True)
                    with open(dest,'wb') as out: shutil.copyfileobj(t.extractfile(m),out,8*1024*1024)
                    print('Extracted',m.name,m.size,flush=True)
        while b:=gz.read(8*1024*1024): total+=len(b)
    (AUDIT/(archive+'.validation.json')).write_text(json.dumps({'bytes':path.stat().st_size,'sha256':h.hexdigest(),'tar_readable':True,'gzip_crc_pass':True,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2))
save_rows(AUDIT/'scct_archive_members.csv',inventory)
# CRC RDS remote size vs downloaded bytes. HEAD failures remain unknown, not complete.
manifest=list(csv.DictReader(open(SOURCE/'CellResDB'/'CellResDB_download_manifest.csv',encoding='utf-8-sig')))
def head(row):
    p=SOURCE/'CellResDB'/row['file']; d={'file':row['file'],'local_bytes':p.stat().st_size if p.exists() else 0,'url':row['url']}
    try:
        with urllib.request.urlopen(urllib.request.Request(row['url'],method='HEAD'),timeout=60) as r:
            d['remote_bytes']=r.headers.get('Content-Length',''); d['etag']=r.headers.get('ETag','');d['status']='size_match' if str(d['local_bytes'])==d['remote_bytes'] else 'size_unknown_or_mismatch'
    except Exception as e: d.update(remote_bytes='',etag='',status=str(e))
    return d
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    rr=list(pool.map(head,[r for r in manifest if '_CRC_' in r['file'] or '_Colorectal ' in r['file']]))
save_rows(AUDIT/'cellres_remote_sizes.csv',rr)
print('Archive audit completed',flush=True)
