"""Extract official CellResDB sample tables and GEO original cell annotations.
Patient response is kept by endpoint; inconsistent labels are never collapsed.
"""
from html.parser import HTMLParser
import pathlib,csv,re,json,gzip,shlex
ROOT=pathlib.Path(__file__).resolve().parents[1]
class Tables(HTMLParser):
    def __init__(self):super().__init__();self.tables=[];self.t=None;self.row=None;self.cell=None
    def handle_starttag(self,tag,attrs):
        if tag=='table':self.t=[]
        if tag=='tr' and self.t is not None:self.row=[]
        if tag in ('td','th') and self.row is not None:self.cell=[]
    def handle_data(self,s):
        if self.cell is not None:self.cell.append(s)
    def handle_endtag(self,tag):
        if tag in ('td','th') and self.cell is not None:self.row.append(' '.join(' '.join(self.cell).split()));self.cell=None
        if tag=='tr' and self.row is not None:self.t.append(self.row);self.row=None
        if tag=='table' and self.t is not None:self.tables.append(self.t);self.t=None
records=[]
for f in sorted((ROOT/'audit/official').glob('*.rds.html')):
    p=Tables();p.feed(f.read_text(encoding='utf-8',errors='replace'))
    sample=[t for t in p.tables if t and 'Patient ID' in t[0]]
    if not sample:continue
    for row in sample[0][1:]:
        if len(row)<5:continue
        r=dict(zip(sample[0][0],row));r['dataset']=f.name[:-9];records.append(r)
fields=list(dict.fromkeys(k for r in records for k in r))
with open(ROOT/'audit/official_sample_tables.csv','w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)
print(fields);print(records[:2]);print('Sample table rows',len(records))
