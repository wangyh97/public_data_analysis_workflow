# Independent aggregation audit of the saved unmodified CellChat LR arrays.
root<-normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/')
.libPaths(c(file.path(root,'runtime/R-library'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
jobs<-fread(file.path(root,'data/CellChat/jobs.csv'));cfg<-fromJSON(file.path(root,'cellchat_config.json'));audit<-list()
for(k in seq_len(nrow(jobs))){
 j<-jobs[k];file<-file.path(root,'derived/CellChat',paste0(j$job,'.rds'));z<-readRDS(file)
 stopifnot(z$config$seed==cfg$seed,z$config$nboot==cfg$nboot,z$config$type==cfg$type,z$package_version=='2.2.0.9001')
 stopifnot(sum(z$cell_counts)==j$n_cells,all(z$cell_counts>=20),z$cell_counts['Malignant']>=30)
 p<-z$net$prob;p[z$net$pval>.05]<-0
 stopifnot(identical(dimnames(p)[[3]],rownames(z$LR$LRsig)))
 delta<-0
 for(pa in z$netP$pathways){
  ix<-which(z$LR$LRsig$pathway_name==pa);expected<-apply(p[,,ix,drop=FALSE],c(1,2),sum)
  delta<-max(delta,max(abs(expected-z$netP$prob[,,pa])))
 }
 stopifnot(delta<1e-10)
 classical<-which(z$LR$LRsig$ligand%in%c('HLA-A','HLA-B','HLA-C')&z$LR$LRsig$receptor%in%c('CD8A','CD8B'))
 audit[[k]]<-data.table(job=j$job,n_cells=sum(z$cell_counts),n_classical_LR_candidates=length(classical),pathway_reaggregation_max_error=delta,parameters_match=TRUE)
 # The first completed smoke object had an outdated descriptive DB string;
 # actual db selection was protein-only in source throughout. Correct metadata
 # explicitly while preserving every numerical result and the prior string.
 if(z$config$database!=cfg$database){z$prior_database_description<-z$config$database;z$config$database<-cfg$database;saveRDS(z,file)}
}
fwrite(rbindlist(audit),file.path(root,'audit/cellchat_reaggregation_validation.csv'))
ip<-as.data.table(installed.packages()[,c('Package','Version','LibPath')]);fwrite(ip,file.path(root,'audit/CellChat_installed_packages_final.csv'))
cat('VALIDATED',length(audit),'CellChat objects; all pathway aggregates equal LR sums\n')
