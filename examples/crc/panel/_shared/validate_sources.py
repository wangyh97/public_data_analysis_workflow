"""Independent non-visual verification: original extraction vs re-extraction,
patient-level rho/BH calculations, source counts, PDF geometry/text, and hashes.
Uses Python solely for numerical/file inspection, never figure generation.
Run: python validate_sources.py
"""
import pathlib,json,hashlib,math
import pandas as pd
import numpy as np
import pdfplumber
ROOT=pathlib.Path(__file__).resolve().parents[2]
PAN=ROOT/'panel';result=[]
ids=['GEO_GSE236581','scCT_dataset_123-I','scCT_dataset_124-I','scCT_dataset_125-I']
for ds in ids:
    old=ROOT/'derived'/f'{ds}_expression.csv.gz'
    new=ROOT/'derived/panel_v2'/f'{ds}_expression.csv.gz'
    cols=list(pd.read_csv(new,nrows=0).columns)
    shared=[c for c in pd.read_csv(old,nrows=0).columns if c in cols]
    readers=[pd.read_csv(f,usecols=shared,chunksize=50000) for f in [old,new]]
    n=0;maxdiff=0.
    for a,b in zip(*readers,strict=True):
        assert list(a.barcode)==list(b.barcode)
        num=[c for c in shared if c!='barcode']
        delta=np.abs(a[num].to_numpy()-b[num].to_numpy()).max()
        maxdiff=max(maxdiff,float(delta));n+=len(a)
    assert maxdiff<1e-9
    result.append(dict(check='full_reextraction_matches_previous',dataset=ds,cells=n,columns=len(shared),max_difference=maxdiff))
for panel in [PAN/'01_MHC_I_correlation',PAN/'02_T_cell_state_correlation']:
    pts=pd.read_csv(panel/'source/patient_values.csv');st=pd.read_csv(panel/'source/correlation_statistics.csv')
    assert not pts.duplicated(['cohort','patient','target','outcome']).any()
    for row in st.itertuples():
        d=pts[(pts.cohort==row.cohort)&(pts.target==row.target)&(pts.outcome==row.outcome)]
        rho=np.corrcoef(d.x.rank(),d.y.rank())[0,1]
        assert len(d)==row.n and abs(rho-row.rho)<1e-12
    for _,d in st.groupby(['cohort','family']):
        order=np.argsort(d.p);p=d.p.to_numpy()[order]
        q=np.minimum(1,np.minimum.accumulate((p*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
        assert np.max(np.abs(q-d.q.to_numpy()[order]))<1e-12
    result.append(dict(check='independent_rho_and_BH',panel=panel.name,tests=len(st),patient_pairs=len(pts)))
    for f in sorted((panel/'figure').glob('*.pdf')):
        with pdfplumber.open(f) as doc:
            assert len(doc.pages)==1
            p=doc.pages[0];chars=p.chars
            assert chars
            outside=[c for c in chars if c['x0']<-.5 or c['x1']>p.width+.5 or c['top']<-.5 or c['bottom']>p.height+.5]
            assert not outside,(f.name,outside[:1])
            # pdfplumber's `size` is a vertical bounding extent for rotated
            # characters, NOT the font size (e.g. a 9 pt rotated 'i' reports 2).
            # Cairo encodes its font scale in the text matrix. Its column norm
            # is rotation invariant and measures the effective glyph em size.
            floor=min(math.hypot(c['matrix'][2],c['matrix'][3]) for c in chars)
            assert floor>=4.99,(f.name,floor)
            result.append(dict(check='PDF_text_geometry',file=str(f.relative_to(PAN)),minimum_font_pt=round(floor,3),outside_page=0))
    # Explicit missingness/exclusion and patient eligibility remain in CSVs.
    eligible=pd.read_csv(panel/'source/sample_compartment_expression.csv')
    result.append(dict(check='cohort_sizes',panel=panel.name,patients=pts.groupby('cohort').patient.nunique().to_dict(),cells_by_compartment=eligible.groupby(['cohort','compartment']).n_cells.sum().rename('cells').reset_index().to_dict('records')))
(PAN/'QA_numeric_and_PDF.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
inputs=list((ROOT/'data/GEO').glob('GSE236581*'))
inputs=[f for f in inputs if 'VDJ' not in f.name]
for ds in ids[1:]:inputs.extend((ROOT/'data/scCT-DB'/ds).glob('*'))
inputs.extend([ROOT/'audit/official_sample_tables.csv',ROOT/'audit/scct_crc_samples.csv',PAN/'_shared/extraction_genes.json'])
manifest=[]
for f in inputs:
    digest=hashlib.sha256()
    with f.open('rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b''):digest.update(chunk)
    manifest.append(dict(path=str(f.relative_to(ROOT)),bytes=f.stat().st_size,sha256=digest.hexdigest()))
for panel in [PAN/'01_MHC_I_correlation',PAN/'02_T_cell_state_correlation']:pd.DataFrame(manifest).to_csv(panel/'source/input_manifest.csv',index=False)
print('Input manifest hashes completed')
