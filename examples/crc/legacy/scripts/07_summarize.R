# Per-sample compartment summaries. Independent unit = patient, never individual cell.
# Keeping portal-specific results allows sensitivity comparisons without counting
# the same GEO study as an independent validation cohort.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
root<-normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/')
cfg<-fromJSON(file.path(root,'analysis_config.json'));set.seed(cfg$seed)
dir.create(file.path(root,'tables'),showWarnings=FALSE);dir.create(file.path(root,'figures'),showWarnings=FALSE)
meta_sc<-fread(file.path(root,'audit/scct_crc_samples.csv'))
sets<-cfg$gene_sets;targets<-cfg$targets
classify<-function(x){
  x<-tolower(x);out<-rep('Other',length(x))
  out[grepl('fibro|stromal',x)]<-'Fibroblast';out[grepl('endothel',x)]<-'Endothelial'
  out[grepl('myelo|monoc|macro|dendritic|neutro|mast',x)]<-'Myeloid'
  out[grepl('b cell|b_cell|plasma',x)]<-'B_plasma'
  out[grepl('nk|natural killer|innate lymphocyte|^ilc$',x)]<-'NK_ILC'
  out[grepl('t cell|t_cell|cd8|cd4|treg|tex|ttr|trm',x)]<-'T'
  out[grepl('epithel',x)]<-'Epithelial_uncertain'
  out[grepl('malignant|cancer cell|tumor cell',x)&!grepl('non.malignant',x)]<-'Malignant'
  out
}
files<-list.files(file.path(root,'derived'),pattern='_expression.csv.gz$',full.names=TRUE)
# Native R checkpoints for small CellRes datasets remain canonical and avoid
# duplicate use of stream validation files.
files<-files[!grepl('GSE205506',files)]
files<-c(list.files(file.path(root,'derived'),pattern='GSE205506.*\\.rds$',full.names=TRUE),files)
summaries<-list();coverage<-list();within<-list();qcs<-list();sample_records<-list();seen<-character()
for(f in files){
 id<-sub('_expression.csv.gz$|.rds$','',basename(f));cat('SUMMARIZE',id,'\n');flush.console()
 if(endsWith(f,'.rds')){
   z<-readRDS(f);m<-as.data.table(z$metadata,keep.rownames='barcode');e<-as.data.table(t(z$counts));e[,barcode:=m$barcode];e[,library:=z$library];e[,nFeature:=m$nFeature_RNA];e[,percent_mt:=m$percent.mt];rm(z)
 }else{e<-fread(f);m<-if(startsWith(id,'scCT'))fread(file.path(root,'data/scCT-DB',id,'Cell info.txt')) else fread(file.path(root,'derived',paste0(id,'_metadata.csv.gz')))}
 if(startsWith(id,'scCT')){
   setnames(m,c('Barcode','Cell type','Sample ID'),c('barcode','celltype','sample'))
   s<-unique(meta_sc[`scCT-DB dataset ID`==id],by='Sample ID');idx<-match(m$sample,s[['Sample ID']]);stopifnot(!anyNA(idx))
   m[,`:=`(portal='scCT-DB',study=s[['Source dataset']][idx],patient=s[['Patient ID']][idx],time=ifelse(s[['Sampling time']][idx]=='Baseline','Pre','Post'),location=s[['Sampling location category']][idx],response=ifelse(s[['Response group category']][idx]=='Response','R','NR'),regimen=s[['Therapeutic regimen']][idx])]
 }else{
   m[,`:=`(portal=if(startsWith(id,'GEO_'))'GEO' else 'CellResDB',study=if(startsWith(id,'GEO_'))sub('GEO_','',id) else sub('_.*','',id),sample=as.character(orig.ident),celltype=as.character(cluster),time=as.character(treatment),location=as.character(Sample_Location),response=as.character(Response),regimen=as.character(drug_name))]
   if('Ident'%in%names(m)){
     m[,sample:=as.character(Ident)]
     m[,time:=ifelse(Treatment=='I','Pre',paste0('Post_',Treatment))]
     m[SubCellType=='c91_Epi_Tumor',celltype:='Malignant epithelial cell']
   }
   m[,patient:=paste(study,patient,sep='-')]
 }
 m[,dataset:=id];e<-e[match(m$barcode,e$barcode)];stopifnot(identical(e$barcode,m$barcode))
 if('nCount_RNA'%in%names(m)){
   dif<-e$library-as.numeric(m$nCount_RNA)
   fwrite(data.table(dataset=id,n_cells=length(dif),max_abs_library_difference=max(abs(dif),na.rm=TRUE),fraction_library_equal=mean(abs(dif)<1e-8,na.rm=TRUE)),file.path(root,'audit',paste0(id,'_library_vs_metadata.csv')))
 }
 # Duplicated samples within a portal/study are removed as units, not mixed.
 keys<-paste(m$portal,m$study,m$sample,sep='|');keep<-!keys%in%seen;ndup<-sum(!keep);seen<-union(seen,unique(keys))
 ncancer<-0L
 if('is_crc'%in%names(m)){ncancer<-sum(keep & !m$is_crc);keep<-keep & m$is_crc}
 coverage[[id]]<-data.table(dataset=id,n_input=nrow(m),n_excluded_duplicate_sample=ndup,n_excluded_nonCRC_or_unmapped=ncancer,n_kept=sum(keep),n_genes_target_available=sum(targets%in%names(e)))
 m<-m[keep];e<-e[keep];if(!nrow(m))next
 m[,compartment:=classify(celltype)]
 m[,tumor_tissue:=tolower(location)%in%c('tumor','colon','rectum','primary tumor','primary tumour')]
 sample_records[[id]]<-unique(m[,.(portal,study,dataset,sample,patient,time,location,response,regimen,tumor_tissue)])
 g<-intersect(c(targets,unique(unlist(sets)), 'CD8A','CD8B','CD3D','CD274','PDCD1LG2','PVR','NECTIN2','LGALS9','CXCL9','CXCL10','CXCR3','CD80','CD86','CD28'),names(e))
 counts<-as.matrix(e[,..g]);norm<-log1p(counts/pmax(e$library,1)*10000)
 sc<-sapply(sets,function(gs)rowMeans(norm[,intersect(gs,g),drop=FALSE]));colnames(sc)<-names(sets)
 group<-interaction(m$sample,m$compartment,drop=TRUE)
 groups<-split(seq_len(nrow(m)),group)
 cd8<-if('SubCellType'%in%names(m))grepl('CD8',m$SubCellType) else grepl('CD8',m$celltype,ignore.case=TRUE)
 if(any(cd8)){
   extra<-split(which(cd8),m$sample[cd8]);names(extra)<-paste0('CD8::',names(extra));groups<-c(groups,extra)
 }
 for(gn in names(groups)){
   ix<-groups[[gn]]
   row<-m[ix[1],.(portal,study,dataset,sample,patient,time,location,response,regimen,tumor_tissue,compartment)]
   if(startsWith(gn,'CD8::'))row[,compartment:='CD8_T']
   row[,`:=`(n_cells=length(ix),total_library=sum(e$library[ix]),mean_nFeature=mean(e$nFeature[ix]),mean_percent_mt=mean(e$percent_mt[ix]))]
   for(j in seq_along(g)){
     row[[paste0('mean_',g[j])]]<-mean(norm[ix,j]);row[[paste0('pb_',g[j])]]<-log2(1+sum(counts[ix,j])/max(sum(e$library[ix]),1)*1e6);row[[paste0('pct_',g[j])]]<-mean(counts[ix,j]>0)
   }
   for(s in names(sets))row[[paste0('score_',s)]]<-mean(sc[ix,s])
   summaries[[length(summaries)+1]]<-row
   if(row$compartment%in%c('Malignant','Epithelial_uncertain')&&length(ix)>=cfg$min_tumor_cells){
     for(tg in targets)for(s in c('MHC_I_core','APM')){
       xx<-norm[ix,tg]; yy<-sc[ix,s]
       if(sd(xx)>0&&sd(yy)>0)within[[length(within)+1]]<-data.table(portal=row$portal,study=row$study,sample=row$sample,patient=row$patient,time=row$time,tumor_tissue=row$tumor_tissue,compartment=row$compartment,gene=tg,signature=s,n_cells=length(ix),rho=cor(xx,yy,method='spearman'))
     }
   }
 }
 qcs[[id]]<-data.table(dataset=id,metric=rep(c('library','nFeature','percent_mt'),each=7),quantile=rep(c(0,.01,.1,.5,.9,.99,1),3),value=c(quantile(e$library,c(0,.01,.1,.5,.9,.99,1)),quantile(e$nFeature,c(0,.01,.1,.5,.9,.99,1)),quantile(e$percent_mt,c(0,.01,.1,.5,.9,.99,1))))
 rm(e,m,counts,norm,sc);gc()
}
a<-rbindlist(summaries,fill=TRUE);samples<-unique(rbindlist(sample_records,fill=TRUE))
# Propagate the final observed R/NR and post-treatment regimen only within patient
# and portal. Conflicting endpoints are explicitly left unresolved.
mapping<-samples[,.(final_response={v<-unique(response[response%in%c('R','NR')]);if(length(v)==1)v else NA_character_},response_conflict=uniqueN(response[response%in%c('R','NR')])>1,
                  final_regimen={v<-unique(regimen[grepl('^Post',time)&!regimen%in%c('None','NA','')]);if(length(v)==1)v else paste(sort(v),collapse=' | ')}),by=.(portal,study,patient)]
a<-merge(a,mapping,by=c('portal','study','patient'),all.x=TRUE);samples<-merge(samples,mapping,by=c('portal','study','patient'),all.x=TRUE)
a[,response_filled:=!response%in%c('R','NR')&!is.na(final_response)]
fwrite(a,file.path(root,'tables/sample_compartment_summary.csv'));fwrite(samples,file.path(root,'tables/sample_manifest.csv'));fwrite(mapping,file.path(root,'tables/patient_response_mapping.csv'))
fwrite(rbindlist(coverage),file.path(root,'audit/analysis_cell_coverage.csv'));fwrite(rbindlist(qcs),file.path(root,'tables/qc_quantiles.csv'));fwrite(rbindlist(within),file.path(root,'tables/within_sample_correlations.csv'))
writeLines(capture.output(sessionInfo()),file.path(root,'audit/R_analysis_sessionInfo.txt'))
