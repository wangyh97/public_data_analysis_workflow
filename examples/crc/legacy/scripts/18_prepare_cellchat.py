"""Freeze CellChat sample eligibility and deterministic stratified cell selection.
Full gene matrices are subsequently extracted; high/low labels are NEVER used
for selecting cells or defining cell types. Author QC and labels are reused.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/CellChat';OUT.mkdir(exist_ok=True)
cfg=dict(seed=20260905,max_cells_per_sample_group=300,min_group_cells=20,
         min_malignant_cells=30,nboot=100,type='triMean',population_size=False,
         raw_use=True,cohort_stage='Pre',database='CellChatDB.human; protein signaling only (exclude Non-protein Signaling)',
         inference_unit='original tumor sample; equal sample mean within patient',
         test_unit='patient',min_patients_each_high_low=3,permutations=10000)
(ROOT/'cellchat_config.json').write_text(json.dumps(cfg,indent=2))
groups=pd.read_csv(ROOT/'tables/patient_high_low_groups.csv')
def classify(label):
 s=label.lower()
 if 'malignant' in s or 'c91_epi_tumor' in s:return 'Malignant'
 if 'treg' in s:return 'Treg'
 if 'cd8' in s:return 'CD8_T'
 if 'cd4' in s:return 'CD4_T'
 if any(k in s for k in ['_t_mki','gdt','dnt','dpt']):return 'Other_T'
 if any(k in s for k in ['mono','mph','macro','mast','neu_','dc_']):return 'Myeloid'
 if any(k in s for k in ['naiveb','memb','gcb','plasmab','b cell','plasma cell']):return 'B_plasma'
 if any(k in s for k in ['endo_','endothelial']):return 'Endothelial'
 if any(k in s for k in ['fibro','pericyte','smc_']):return 'Fibroblast'
 if 'ilc' in s or 'nk ' in s:return 'NK_ILC'
 return 'Other'
m=pd.read_csv(ROOT/'derived/GEO_GSE236581_metadata.csv.gz')
m['matrix_col']=np.arange(len(m))+1
m['patient']='GSE236581-'+m.patient
m=m[(m.Treatment=='I')&(m.Tissue=='Tumor')&m.is_crc&m.patient.isin(groups.patient)].copy()
m['label']=m.SubCellType.map(classify);m['sample']=m.Ident;m['study']='GSE236581';m['dataset']='GEO'
allmeta=[m[['barcode','matrix_col','patient','label','sample','study','dataset']]]
sm=pd.read_csv(ROOT/'audit/scct_crc_samples.csv')
for ds in ['scCT_dataset_123-I','scCT_dataset_124-I','scCT_dataset_125-I']:
 c=pd.read_csv(ROOT/'data/scCT-DB'/ds/'Cell info.txt',sep='\t').rename(columns={'Barcode':'barcode','Cell type':'celltype','Sample ID':'sample'})
 s=sm[(sm['scCT-DB dataset ID']==ds)&(sm['Sampling time']=='Baseline')]
 c=c[c['sample'].isin(s['Sample ID'])].copy();c['patient']=c['sample'].map(s.set_index('Sample ID')['Patient ID'])
 c=c[c.patient.isin(groups.patient)];c['label']=c.celltype.map(classify)
 c['study']='GSE205506';c['dataset']=ds;c['matrix_col']=0
 allmeta.append(c[['barcode','matrix_col','patient','label','sample','study','dataset']])
m=pd.concat(allmeta,ignore_index=True)
sizes=m.groupby(['study','dataset','sample','patient','label']).size().rename('n_available').reset_index()
sizes['eligible_group']=sizes.n_available>=np.where(sizes.label=='Malignant',30,20)
ok=sizes[(sizes.label=='Malignant')&sizes.eligible_group]['sample']
m=m[m['sample'].isin(ok)&(m.label!='Other')]
rng=np.random.default_rng(cfg['seed']);chosen=[]
for key,g in m.groupby(['study','dataset','sample','patient','label'],sort=True):
 if len(g)<(30 if key[-1]=='Malignant' else 20):continue
 chosen.append(g.iloc[np.sort(rng.choice(len(g),min(len(g),300),replace=False))])
sel=pd.concat(chosen).sort_values(['study','dataset','sample','label','barcode']).reset_index(drop=True)
sel['job']=sel.study+'__'+sel['sample'];sel['local_col']=sel.groupby('job').cumcount()+1
sel.to_csv(OUT/'selected_cells.csv',index=False)
sizes.to_csv(ROOT/'audit/cellchat_cell_group_eligibility.csv',index=False)
jobs=sel.groupby(['job','study','dataset','sample','patient']).size().rename('n_cells').reset_index()
jobs.to_csv(OUT/'jobs.csv',index=False)
for job,g in sel.groupby('job'):g.to_csv(OUT/(job+'_cells.csv'),index=False)
print(jobs.to_string(index=False));print('TOTAL',len(sel),'cells',len(jobs),'samples')
