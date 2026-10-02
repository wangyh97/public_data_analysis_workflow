# Shared processing engine. Both panel-specific data_processing.R entry points
# call this file. Only deposited processed counts and published annotations are
# used; no new CNV, doublet, ambient-RNA or hard-QC calls are claimed.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
process_panel <- function(panel_dir, rebuild_raw=FALSE) {
 root <- normalizePath(file.path(panel_dir,'../..'),winslash='/')
 cfg <- fromJSON(file.path(panel_dir,'script/config.json'))
 paths <- fromJSON(file.path(root,'analysis_config.json'))$paths
 cache <- file.path(root,'derived/panel_v2');dir.create(cache,recursive=TRUE,showWarnings=FALSE)
 # Node performs a bounded-memory integer Matrix Market read; Python streams the
 # gene-by-cell scCT CSV. Every gene contributes to each cell's library size.
 if(rebuild_raw || !file.exists(file.path(cache,'GEO_GSE236581_validation.json'))) {
   status <- system2(paths$node,c(shQuote(file.path(root,'panel/_shared/extract_geo.mjs')),shQuote(root)))
   stopifnot(status==0)
 }
 if(rebuild_raw) for(id in c('123','124','125')) unlink(file.path(cache,paste0('scCT_dataset_',id,'-I.json')))
 status <- system2(paths$python,c(shQuote(file.path(root,'panel/_shared/extract_scct.py')),shQuote(root)))
 stopifnot(status==0)
 all_genes<-fromJSON(file.path(root,'panel/_shared/extraction_genes.json'))
 stopifnot(all(unique(c(cfg$targets,unlist(cfg$gene_sets))) %in% all_genes))
 set.seed(cfg$seed);summaries<-list();audit<-list();coverage<-list()
 # CRC membership is frozen from the official CRC sample table; no clinical
 # response is required for these correlations and R/NR is never pooled.
 clinical<-fread(file.path(root,'audit/official_sample_tables.csv'))
 crc_ids<-unique(clinical[grepl('GSE236581',dataset)][['Patient ID']])
 smeta<-fread(file.path(root,'audit/scct_crc_samples.csv'))
 ids<-c('GEO_GSE236581','scCT_dataset_123-I','scCT_dataset_124-I','scCT_dataset_125-I')
 seen<-character()
 for(id in ids) {
   message('Processing ',id)
   if(startsWith(id,'GEO')) {
     m<-as.data.table(read.table(gzfile(file.path(root,'data/GEO/GSE236581_CRC-ICB_metadata.txt.gz')),header=TRUE,row.names=1),keep.rownames='barcode')
     m[,`:=`(cohort='GSE236581',patient=paste0('GSE236581-',Patient),sample=Ident,
       annotation=SubCellType,malignant=SubCellType=='c91_Epi_Tumor',
       T_cell=MajorCellType=='T',CD8=MajorCellType=='T' & grepl('CD8',SubCellType),
       eligible_tissue=tolower(Tissue)%in%c('tumor','colon','rectum','primary tumor'),
       eligible_time=Treatment=='I',eligible_crc=Patient%in%crc_ids,
       annotation_source='GEO author c91_Epi_Tumor')]
   } else {
     m<-fread(file.path(root,'data/scCT-DB',id,'Cell info.txt'))
     setnames(m,c('Barcode','Cell type','Sample ID'),c('barcode','annotation','sample'))
     s<-unique(smeta[get('scCT-DB dataset ID')==id],by='Sample ID')
     ix<-match(m$sample,s[['Sample ID']]);stopifnot(!anyNA(ix))
     m[,`:=`(cohort='GSE205506',patient=s[['Patient ID']][ix],
       malignant=annotation=='Malignant epithelial cell',
       T_cell=grepl('t cell|t_cell|cd8|cd4|treg|tex|ttr|trm',annotation,ignore.case=TRUE),
       CD8=grepl('CD8',annotation,ignore.case=TRUE),
       eligible_tissue=tolower(s[['Sampling location category']][ix])=='tumor',
       eligible_time=s[['Sampling time']][ix]=='Baseline',eligible_crc=TRUE,
       annotation_source='scCT-DB malignant label; not independently CNV verified')]
   }
   keys<-paste(m$cohort,m$sample);m[,duplicate_sample:=keys%in%seen];seen<-union(seen,unique(keys))
   audit[[id]]<-m[,.(n_cells=.N),by=.(cohort,patient,sample,eligible_tissue,eligible_time,eligible_crc,duplicate_sample)]
   m<-m[eligible_tissue & eligible_time & eligible_crc & !duplicate_sample]
   e<-fread(file.path(cache,paste0(id,'_expression.csv.gz')))
   missing<-setdiff(all_genes,names(e))
   coverage[[id]]<-data.table(dataset=id,gene=all_genes,present=all_genes%in%names(e))
   if(length(missing))stop('Missing requested genes in ',id,': ',paste(missing,collapse=','))
   stopifnot(!anyDuplicated(e$barcode),!anyDuplicated(m$barcode))
   e<-e[match(m$barcode,e$barcode)];stopifnot(identical(m$barcode,e$barcode),all(e$library>0))
   norm<-log1p(as.matrix(e[,..all_genes])/e$library*10000)
   for(comp in c('Malignant','All_T','CD8_T')) {
     ix<-which(switch(comp,Malignant=m$malignant,All_T=m$T_cell,CD8_T=m$CD8))
     for(j in split(ix,m$sample[ix])) {
       z<-m[j[1],.(cohort,patient,sample,annotation_source)]
       z[,`:=`(compartment=comp,n_cells=length(j))]
       for(g in all_genes)z[[g]]<-mean(norm[j,g])
       summaries[[length(summaries)+1]]<-z
     }
   }
   rm(e,m,norm);gc()
 }
 a<-rbindlist(summaries);src<-file.path(panel_dir,'source')
 fwrite(rbindlist(audit),file.path(src,'input_sample_audit.csv'))
 fwrite(rbindlist(coverage),file.path(src,'gene_coverage.csv'))
 # Scores are arithmetic means of listed genes' mean log1pCP10K; by linearity
 # this equals mean per-cell module expression. No control genes or z-scoring.
 for(s in names(cfg$gene_sets))set(a,j=paste0('score_',s),value=rowMeans(as.matrix(a[,cfg$gene_sets[[s]],with=FALSE])))
 a[,eligible_compartment:=n_cells>=ifelse(compartment=='Malignant',cfg$minimum_tumor_cells,cfg$minimum_T_cells)]
 fwrite(a,file.path(src,'sample_compartment_expression.csv'))
 definitions<-rbindlist(lapply(names(cfg$gene_sets),function(s)data.table(gene_set=s,gene=cfg$gene_sets[[s]])))
 fwrite(definitions,file.path(src,'gene_sets.csv'))
 if(cfg$panel_type=='MHC') {
   spec<-rbind(data.table(outcome=names(cfg$gene_sets),column=paste0('score_',names(cfg$gene_sets)),compartment='Malignant',family='MHC_modules'),
     data.table(outcome=cfg$gene_sets$MHC_I_structure,column=cfg$gene_sets$MHC_I_structure,compartment='Malignant',family='MHC_structural_single_genes'))
 } else {
   spec<-rbind(data.table(outcome=paste0('All_T_',names(cfg$gene_sets)),column=paste0('score_',names(cfg$gene_sets)),compartment='All_T',family='T_cell_states'),
     data.table(outcome=paste0('CD8_T_',names(cfg$gene_sets)[1:3]),column=paste0('score_',names(cfg$gene_sets)[1:3]),compartment='CD8_T',family='T_cell_states'))
 }
 records<-list()
 for(k in seq_len(nrow(spec)))for(g in cfg$targets) {
   ss<-spec[k]
   x<-a[compartment=='Malignant' & eligible_compartment,.(cohort,patient,sample,x=get(g),n_tumor_cells=n_cells)]
   y<-a[compartment==ss$compartment & eligible_compartment,.(cohort,patient,sample,y=get(ss$column),n_state_cells=n_cells)]
   # Join within the SAME tumor specimen first. Then give eligible specimens
   # equal weight within each patient. This avoids matching unmatched biopsies.
   pair<-merge(x,y,by=c('cohort','patient','sample'))
   pp<-pair[,.(x=mean(x),y=mean(y),n_samples=.N,samples=paste(sort(sample),collapse=';'),n_tumor_cells=sum(n_tumor_cells),n_state_cells=sum(n_state_cells)),by=.(cohort,patient)]
   pp[,`:=`(target=g,outcome=ss$outcome,state_compartment=ss$compartment,family=ss$family)]
   records[[length(records)+1]]<-pp
 }
 pts<-rbindlist(records)
 stopifnot(!anyDuplicated(pts[,.(cohort,patient,target,outcome)]))
 fwrite(pts,file.path(src,'patient_values.csv'))
 stat<-pts[,{
   n<-length(x);valid<-n>=cfg$minimum_patients && sd(x)>0 && sd(y)>0
   if(valid) {
     exact<-n<10 && !anyDuplicated(x) && !anyDuplicated(y)
     ct<-cor.test(x,y,method='spearman',exact=exact)
     boot<-replicate(cfg$bootstrap_replicates,{ii<-sample.int(n,replace=TRUE);if(sd(x[ii])==0||sd(y[ii])==0)NA_real_ else cor(x[ii],y[ii],method='spearman')})
     ci<-quantile(boot,c(.025,.975),na.rm=TRUE,names=FALSE)
     list(n=n,rho=unname(ct$estimate),p=ct$p.value,ci_low=ci[1],ci_high=ci[2],p_method=if(exact)'exact Spearman' else 'asymptotic Spearman',status='tested')
   } else list(n=n,rho=NA_real_,p=NA_real_,ci_low=NA_real_,ci_high=NA_real_,p_method='not tested',status='insufficient patients or constant values')
 },by=.(cohort,target,outcome,state_compartment,family)]
 stat[,q:=p.adjust(p,method='BH'),by=.(cohort,family)]
 fwrite(stat,file.path(src,'correlation_statistics.csv'))
 fwrite(spec,file.path(src,'outcome_definitions.csv'))
 writeLines(capture.output(sessionInfo()),file.path(src,'R_processing_sessionInfo.txt'))
 message('Completed: ',panel_dir)
}
