# Fig3: raw integer counts -> PATIENT x STAGE pseudobulk -> longitudinal data.
# This is NOT the mean log-normalized expression used in Fig1/2.
# Run Rscript data_processing.R [--rebuild-raw]. Paths are resolved from this file.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
script<-sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])
panel<-normalizePath(file.path(dirname(script),'..'),winslash='/')
root<-normalizePath(file.path(panel,'../..'),winslash='/');src<-file.path(panel,'source')
cfg<-fromJSON(file.path(panel,'script/config.json'));paths<-fromJSON(file.path(root,'analysis_config.json'))$paths
cache<-file.path(root,'derived/panel_v2');rebuild<-'--rebuild-raw'%in%commandArgs(TRUE)
# The shared readers stream every raw count and retain the full-gene library.
# No CellChat downsampling, logarithms, or cell-level normalization occurs here.
if(rebuild || !file.exists(file.path(cache,'GEO_GSE236581_validation.json'))){
 stopifnot(system2(paths$node,c(shQuote(file.path(root,'panel/_shared/extract_geo.mjs')),shQuote(root)))==0)
}
if(rebuild)for(id in c('123','124','125'))unlink(file.path(cache,paste0('scCT_dataset_',id,'-I.json')))
stopifnot(system2(paths$python,c(shQuote(file.path(root,'panel/_shared/extract_scct.py')),shQuote(root)))==0)
clinical<-fread(file.path(root,'audit/official_sample_tables.csv'))
cr<-unique(clinical[grepl('GSE236581',dataset),.(`Patient ID`,Response)])
stopifnot(!anyDuplicated(cr[['Patient ID']]))
sm<-fread(file.path(root,'audit/scct_crc_samples.csv'))
ids<-c('GEO_GSE236581','scCT_dataset_123-I','scCT_dataset_124-I','scCT_dataset_125-I')
samples<-list();audit<-list();coverage<-list();seen<-character()
for(id in ids){
 message('Pseudobulk ',id)
 if(startsWith(id,'GEO')){
  m<-as.data.table(read.table(gzfile(file.path(root,'data/GEO/GSE236581_CRC-ICB_metadata.txt.gz')),header=TRUE,row.names=1),keep.rownames='barcode')
  m[,`:=`(cohort='GSE236581',patient=paste0('GSE236581-',Patient),sample=Ident,
    stage=fifelse(Treatment=='I','Pre',Treatment),raw_stage=Treatment,
    response=cr$Response[match(Patient,cr[['Patient ID']])],
    tumor_tissue=tolower(Tissue)%in%c('tumor','colon','rectum','primary tumor'),
    is_crc=Patient%in%cr[['Patient ID']],malignant=SubCellType=='c91_Epi_Tumor',
    annotation_source='GEO author c91_Epi_Tumor')]
 }else{
  m<-fread(file.path(root,'data/scCT-DB',id,'Cell info.txt'))
  setnames(m,c('Barcode','Cell type','Sample ID'),c('barcode','celltype','sample'))
  s<-unique(sm[get('scCT-DB dataset ID')==id],by='Sample ID');ix<-match(m$sample,s[['Sample ID']]);stopifnot(!anyNA(ix))
  raw_time<-s[['Sampling time']][ix];stopifnot(all(raw_time%in%c('Baseline','Post-treatment')))
  m[,`:=`(cohort='GSE205506',patient=s[['Patient ID']][ix],sample=sample,
    stage=ifelse(raw_time=='Baseline','Pre','Post'),raw_stage=raw_time,
    response=ifelse(s[['Response group category']][ix]=='Response','R',ifelse(s[['Response group category']][ix]%in%c('Non_response','Non-response'),'NR',NA_character_)),
    tumor_tissue=tolower(s[['Sampling location category']][ix])=='tumor',is_crc=TRUE,
    malignant=celltype=='Malignant epithelial cell',
    annotation_source='scCT-DB malignant label; not independently CNV verified')]
 }
 m[,dataset:=id];keys<-paste(m$cohort,m$sample);m[,duplicate_sample:=keys%in%seen];seen<-union(seen,unique(keys))
 audit[[id]]<-m[,.(n_cells=.N,n_malignant=sum(malignant)),by=.(cohort,patient,sample,stage,raw_stage,response,tumor_tissue,is_crc,duplicate_sample,dataset,annotation_source)]
 m<-m[tumor_tissue & is_crc & !duplicate_sample]
 e<-fread(file.path(cache,paste0(id,'_expression.csv.gz')),select=c('barcode','library',cfg$targets))
 coverage[[id]]<-data.table(dataset=id,target=cfg$targets,present=cfg$targets%in%names(e))
 stopifnot(!anyDuplicated(e$barcode),!anyDuplicated(m$barcode))
 e<-e[match(m$barcode,e$barcode)];stopifnot(identical(m$barcode,e$barcode))
 # Preserve a sample record even when no malignant cells were recovered.
 for(ix in split(seq_len(nrow(m)),m$sample)){
  ii<-ix[m$malignant[ix]]
  z<-m[ix[1],.(cohort,patient,sample,stage,raw_stage,response,dataset,annotation_source)]
  z[,`:=`(n_all_tissue_cells=length(ix),n_malignant_cells=length(ii),full_gene_library=sum(e$library[ii]))]
  for(g in cfg$targets)set(z,j=paste0('count_',g),value=sum(e[[g]][ii]))
  samples[[length(samples)+1]]<-z
 }
 rm(m,e);gc()
}
a<-rbindlist(samples);fwrite(a,file.path(src,'sample_pseudobulk_counts.csv'))
fwrite(rbindlist(audit),file.path(src,'input_sample_audit.csv'));fwrite(rbindlist(coverage),file.path(src,'gene_coverage.csv'))
stopifnot(all(a$response%in%c('R','NR')))
stopifnot(a[,uniqueN(response),by=.(cohort,patient)][,all(V1==1)])
# TRUE pseudobulk: add counts AND full-gene libraries across all tumor specimens
# in the same patient/stage. Do not average sample log CPM values.
num<-c('n_all_tissue_cells','n_malignant_cells','full_gene_library',paste0('count_',cfg$targets))
b<-a[,c(lapply(.SD,sum),list(n_samples=.N,samples=paste(sort(sample),collapse=';'))),by=.(cohort,patient,stage,response,annotation_source),.SDcols=num]
b[,eligible:=n_malignant_cells>=cfg$minimum_cells_per_patient_timepoint & full_gene_library>0]
b[,exclusion_reason:=fifelse(eligible,'included',fifelse(n_malignant_cells==0,'no malignant-labelled cells','fewer than minimum malignant-labelled cells'))]
for(g in cfg$targets)set(b,j=paste0('CPM_',g),value=ifelse(b$full_gene_library>0,b[[paste0('count_',g)]]/b$full_gene_library*1e6,NA_real_))
long<-rbindlist(lapply(cfg$targets,function(g){z<-copy(b);z[,`:=`(gene=g,raw_count=get(paste0('count_',g)),CPM=get(paste0('CPM_',g)))];z[,expression:=log2(1+CPM)];z[,c(paste0('count_',cfg$targets),paste0('CPM_',cfg$targets)):=NULL];z}))
elig<-unique(long[eligible==TRUE,.(cohort,patient,stage,response)])[,.(has_pre=any(stage=='Pre'),n_post=sum(stage!='Pre')),by=.(cohort,patient,response)]
elig[,trajectory_included:=has_pre & n_post>0]
long<-merge(long,elig,by=c('cohort','patient','response'),all.x=TRUE)
long[is.na(trajectory_included),trajectory_included:=FALSE]
long[,stage_order:=ifelse(cohort=='GSE236581',match(stage,c('Pre','II','III','IV')),match(stage,c('Pre','Post')))]
stopifnot(!anyNA(long$stage_order),!anyDuplicated(long[,.(cohort,patient,stage,gene)]))
setorder(long,cohort,patient,gene,stage_order)
fwrite(long,file.path(src,'patient_pseudobulk_all.csv'))
plotdata<-long[eligible & trajectory_included]
fwrite(plotdata,file.path(src,'trajectory_plot_data.csv'));fwrite(elig,file.path(src,'patient_eligibility.csv'))
availability<-unique(long[,.(cohort,patient,stage,response,eligible,trajectory_included)])[,.(patients_available=.N,patients_eligible=sum(eligible),patients_in_trajectories=sum(eligible & trajectory_included)),by=.(cohort,stage,response)]
fwrite(availability,file.path(src,'timepoint_availability.csv'))
base<-plotdata[stage=='Pre',.(cohort,patient,gene,baseline_expression=expression)]
delta<-merge(plotdata[stage!='Pre'],base,by=c('cohort','patient','gene'))
delta[,delta_log2CPM:=expression-baseline_expression]
fwrite(delta,file.path(src,'change_from_baseline.csv'))
# Exploratory paired tests. R/NR labels are not pooled across studies/endpoints.
tests<-rbindlist(lapply(c('All','R','NR'),function(gr){
 d<-if(gr=='All')delta else delta[response==gr]
 if(!nrow(d))return(NULL)
 d[,{
  valid<-.N>=cfg$minimum_pairs_for_test
  pv<-if(valid && any(delta_log2CPM!=0))suppressWarnings(wilcox.test(delta_log2CPM,mu=0,exact=FALSE)$p.value) else if(valid)1 else NA_real_
  list(group=gr,n_pairs=.N,n_R=sum(response=='R'),n_NR=sum(response=='NR'),median_delta=median(delta_log2CPM),n_increased=sum(delta_log2CPM>0),n_decreased=sum(delta_log2CPM<0),p=pv,status=if(valid)'exploratory paired Wilcoxon' else 'insufficient pairs')
 },by=.(cohort,gene,stage)]
}),fill=TRUE)
tests[,q:=p.adjust(p,'BH'),by=.(cohort,group)]
fwrite(tests,file.path(src,'paired_change_statistics.csv'))
fwrite(unique(a[,.(cohort,stage,raw_stage)]),file.path(src,'timepoint_mapping.csv'))
writeLines(capture.output(sessionInfo()),file.path(src,'R_processing_sessionInfo.txt'))
print(availability);print(tests[group=='All'])
# User-requested primary display statistics. Keep earlier Wilcoxon results as a
# separately named sensitivity table; all new plotted p/q values use paired t.
source(file.path(panel,'script/paired_t_tests.R'))
