"""Read official source pages and HEAD headers; no upload of local data."""
import urllib.request,urllib.parse,pathlib,csv,json,concurrent.futures
ROOT=pathlib.Path(__file__).resolve().parents[1];out=ROOT/'audit'/'official';out.mkdir(exist_ok=True)
SOURCE=pathlib.Path(json.loads((ROOT/'analysis_config.json').read_text())['paths']['source_dir'])
rows=list(csv.DictReader(open(SOURCE/'CellResDB/CellResDB_download_manifest.csv',encoding='utf-8-sig')))
rows=[r for r in rows if '_CRC_' in r['file'] or '_Colorectal ' in r['file']]
def fetch(url,dest,head=False):
    try:
        with urllib.request.urlopen(urllib.request.Request(url,method='HEAD' if head else 'GET',headers={'User-Agent':'CRC-reproducibility-audit/1.0'}),timeout=90) as r:
            h=dict(r.headers)
            if not head:dest.write_bytes(r.read())
            return {'url':url,'headers':h,'status':'ok'}
    except Exception as e:return {'url':url,'status':str(e)}
result=[]
for r in rows:
    h=fetch(r['url'],None,True);h['local_bytes']=(SOURCE/'CellResDB'/r['file']).stat().st_size;h['file']=r['file'];result.append(h)
(out/'cellres_headers.json').write_text(json.dumps(result,indent=2))
urls={**{g+'_GEO.txt':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+g+'&targ=self&form=text&view=full' for g in ['GSE205506','GSE236581','GSE232525']},
      'scct_home.html':'http://scctdb.ncpsb.org.cn/',
      'scct_download.html':'http://scctdb.ncpsb.org.cn/download',
      'cellres_crc.html':'https://cellknowledge.com.cn/cellresponse/php_mysql/getDataSetByOrgan.php?disease=CRC'}
for r in rows:urls[r['file']+'.html']='https://cellknowledge.com.cn/cellresponse/php_mysql/datasetDetail.php?id='+urllib.parse.quote(r['file'][:-4])
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
    rr=list(ex.map(lambda item:fetch(item[1],out/item[0]),urls.items()))
(out/'page_fetch_status.json').write_text(json.dumps(rr,indent=2))
scct_base='http://scctdb.ncpsb.org.cn/biobank/reportTemplate/5e7ac090233167dcfccd9031/group/dashboard/scCT/'
scct=[]
for name in ['SampleInfo.tar.gz','CellInfo.tar.gz','Matrix.tar.gz']:
    h=fetch(scct_base+name,None,True);h['local_bytes']=(SOURCE/'scCT-DB'/name).stat().st_size;scct.append(h)
(out/'scct_archive_headers.json').write_text(json.dumps(scct,indent=2))
print(json.dumps(result,indent=2));print(json.dumps(rr,indent=2))
