"""Independent numerical and PDF-text checks; Python does not draw any figures."""
from pathlib import Path
import pandas as pd
import numpy as np
import pdfplumber,math,json,gzip,csv,hashlib
def paired_t_reference(post,pre):
    # Independent Student-t CDF via Simpson integration after x=sqrt(df)*tan(a).
    # Standard-library math + numpy suffice; no additional installation required.
    delta=np.asarray(post)-np.asarray(pre);n=len(delta);df=n-1
    t=delta.mean()/(delta.std(ddof=1)/np.sqrt(n))
    angle=math.atan(abs(t)/math.sqrt(df));steps=10000
    a=np.linspace(0,angle,steps+1);v=np.cos(a)**(df-1)
    integral=angle/(3*steps)*(v[0]+v[-1]+4*v[1:-1:2].sum()+2*v[2:-1:2].sum())
    norm=math.gamma((df+1)/2)/(math.sqrt(math.pi)*math.gamma(df/2))
    return t,1-2*norm*integral
p=Path(__file__).resolve().parents[1]
a=pd.read_csv(p/'source/sample_pseudobulk_counts.csv')
b=pd.read_csv(p/'source/patient_pseudobulk_all.csv')
d=pd.read_csv(p/'source/trajectory_plot_data.csv')
assert not d.duplicated(['cohort','patient','stage','gene']).any()
for r in b.itertuples():
 s=a[(a.cohort==r.cohort)&(a.patient==r.patient)&(a.stage==r.stage)]
 total=s.full_gene_library.sum();counts=s['count_'+r.gene].sum()
 assert counts==r.raw_count and total==r.full_gene_library
 if total>0:assert abs(np.log2(1+1e6*counts/total)-r.expression)<1e-10
for _,z in d.groupby(['cohort','patient','gene']):
 assert 'Pre' in z.stage.values and len(z)>=2 and (z.n_malignant_cells>=30).all()
assert set(d[d.cohort=='GSE205506'].stage)=={'Pre','Post'}
assert set(d[d.cohort=='GSE236581'].stage)=={'Pre','II','III','IV'}
delta=pd.read_csv(p/'source/change_from_baseline.csv')
assert np.allclose(delta.expression-delta.baseline_expression,delta.delta_log2CPM)
# Original GEO family file confirms that 211 / 213 etc are sample identifiers,
# not longitudinal stage II/III. Match all used GSMs to their treatment labels.
soft=gzip.decompress((p/'source/GSE205506_family.soft.gz').read_bytes()).decode()
catalog=[]
for block in soft.split('^SAMPLE = ')[1:]:
 lines=block.splitlines();r={'sample':lines[0]}
 for line in lines:
  if line.startswith('!Sample_title = '):r['title']=line.split(' = ',1)[1]
  if line.startswith('!Sample_characteristics_ch1 = '):
   key,val=line.split(' = ',1)[1].split(': ',1);r[key]=val
 r['verified_stage']='Pre' if r.get('treatment','').lower()=='untreated' else 'Post'
 catalog.append(r)
cat=pd.DataFrame(catalog);cat.to_csv(p/'source/GSE205506_GEO_sample_catalog.csv',index=False)
used=a[a.cohort=='GSE205506'].merge(cat,on='sample',validate='many_to_one')
assert len(used)==16 and (used.stage==used.verified_stage).all()
assert (used.patient.str.replace('GSE205506-','',regex=False)==used.subject).all()
res={'pseudobulk_rows_checked':len(b),'plotted_gene_points':len(d),
 'patients':d.groupby('cohort').patient.nunique().to_dict(),'GSE205506_GEO_sample_matches':len(used),'figures':[]}
if (p/'source/paired_t_test_statistics.csv').exists():
 ts=pd.read_csv(p/'source/paired_t_test_statistics.csv');pairs=pd.read_csv(p/'source/paired_t_test_pairs.csv')
 checked=0
 for r in ts.itertuples():
  z=pairs[(pairs.cohort==r.cohort)&(pairs.gene==r.gene)&(pairs.stage==r.stage)]
  if r.group!='All':z=z[z.response==r.group]
  assert len(z)==r.n_pairs
  if np.isfinite(r.p):
   t_ref,p_ref=paired_t_reference(z.expression,z.pre);assert abs(p_ref-r.p)<1e-9
   assert abs(t_ref-r.t)<1e-10;checked+=1
 for _,z in ts.groupby(['cohort','group']):
  z=z[z.p.notna()].sort_values('p')
  if len(z):
   adj=np.minimum(1,np.minimum.accumulate((z.p.to_numpy()*len(z)/np.arange(1,len(z)+1))[::-1])[::-1])
   assert np.allclose(z.q,adj,atol=1e-12)
 res['independently_verified_paired_t_tests']=checked
for f in sorted((p/'figure').glob('*.pdf')):
 with pdfplumber.open(f) as doc:
  assert len(doc.pages)==1
  pg=doc.pages[0];chars=pg.chars;assert chars
  outside=[c for c in chars if c['x0']<-.5 or c['x1']>pg.width+.5 or c['top']<-.5 or c['bottom']>pg.height+.5]
  assert not outside,f.name
  # Cairo uses Tf=1 and scales via Tm; account for rotation/scaling.
  size=min(math.hypot(c['matrix'][2],c['matrix'][3]) for c in chars);assert size>=4.99
  res['figures'].append({'file':f.name,'minimum_effective_font_pt':size,'outside_page':0})
(p/'source/QA.json').write_text(json.dumps(res,indent=2))
prior=p.parent/'01_MHC_I_correlation/source/input_manifest.csv'
manifest=pd.read_csv(prior)
for f in [p/'source/GSE205506_family.soft.gz',p/'script/config.json']:
 row=pd.DataFrame([dict(path=str(f.relative_to(p.parents[1])),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest())])
 manifest=pd.concat([manifest,row],ignore_index=True)
manifest.to_csv(p/'source/input_manifest.csv',index=False)
print(json.dumps(res,indent=2))
