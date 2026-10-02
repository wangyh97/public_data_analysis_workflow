# Read each original Seurat RDS without installing/modifying Seurat.
# Only retrieve serialized slots; all matrix math uses installed Matrix.
# Author QC and annotations are preserved. This is a processed-data reanalysis,
# not a fresh raw-barcode QC pipeline; no claim of new doublet removal is made.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(TRUE)
root <- normalizePath(if(length(args)) args[1] else file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'), winslash='/')
src <- file.path(jsonlite::fromJSON(file.path(root,'analysis_config.json'))$paths$source_dir,'CellResDB')
dir.create(file.path(root,'derived'),showWarnings=FALSE)
genes <- c('METTL3','STRAP','PTBP1','HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','TAPBP','NLRC5','PSMB8','PSMB9','PSMB10','ERAP1','ERAP2','CALR','CANX','PDIA3','IRF1','STAT1','ISG15','IFI6','IFIT1','IFIT3','MX1','OAS1','MKI67','TOP2A','PCNA','TYMS','CD3D','CD3E','TRAC','CD8A','CD8B','CD4','NKG7','GNLY','PRF1','GZMA','GZMB','GZMH','IFNG','CCL5','PDCD1','LAG3','HAVCR2','TIGIT','TOX','ENTPD1','CXCL13','TCF7','CCR7','IL7R','LEF1','FOXP3','IL2RA','CTLA4','EPCAM','KRT8','KRT18','KRT19','PTPRC','CD274','PDCD1LG2','PVR','NECTIN2','LGALS9','CXCL9','CXCL10','CXCR3','CD80','CD86','CD28','HLA-DRA','CD74','LST1','CD68','MS4A1')
files <- list.files(src,pattern='GSE205506_CRC_',full.names=TRUE)
if(length(args)>1) files<-files[grepl(args[2],basename(files),fixed=TRUE)]
for (f in files) {
  id <- sub('.rds$','',basename(f)); dest<-file.path(root,'derived',paste0(id,'.rds'))
  if(file.exists(dest)){cat('CHECKPOINT',id,'\n');next}
  cat('READ',id,format(Sys.time()),'\n');flush.console()
  x<-readRDS(f); m<-attr(x,'meta.data'); assay<-attr(x,'assays')[['RNA']]
  counts<-attr(assay,'counts'); stopifnot(inherits(counts,'sparseMatrix'),identical(colnames(counts),rownames(m)))
  stopifnot(!anyDuplicated(rownames(counts)),!anyDuplicated(colnames(counts)),all(is.finite(counts@x)),all(counts@x>=0))
  lib<-Matrix::colSums(counts); present<-intersect(genes,rownames(counts))
  small<-as.matrix(counts[present,,drop=FALSE]); norm<-log1p(sweep(small,2,pmax(lib,1)/10000,'/'))
  # Keep only requested genes per cell. Full matrices remain in source RDS.
  red<-attr(x,'reductions'); embedding<-NULL
  if('umap'%in%names(red)) embedding<-attr(red[['umap']],'cell.embeddings')
  saveRDS(list(id=id,metadata=m,counts=small,lognorm=norm,library=lib,umap=embedding,
               n_genes=nrow(counts),genes_present=present,genes_missing=setdiff(genes,present),
               integer_counts=all(abs(counts@x-round(counts@x))<1e-8)),dest,compress=TRUE)
  keep<-!duplicated(m$sample_ID)
  write.csv(m[keep,,drop=FALSE],file.path(root,'audit',paste0(id,'_samples.csv')),row.names=FALSE)
  write.csv(as.data.frame(table(m$cluster)),file.path(root,'audit',paste0(id,'_celltypes.csv')),row.names=FALSE)
  cat('DONE',id,'genes',nrow(counts),'cells',ncol(counts),'samples',sum(keep),'patients',length(unique(m$patient)),'\n');flush.console()
  rm(x,m,assay,counts,small,norm,red,embedding);gc()
}
writeLines(capture.output(sessionInfo()),file.path(root,'audit','R_sessionInfo.txt'))
