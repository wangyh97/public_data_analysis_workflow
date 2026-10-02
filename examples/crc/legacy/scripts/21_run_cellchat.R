# Real CellChat v2 workflow, one object per original baseline tumor sample.
# The full-gene sparse matrix supplies the normalization denominator. CellChat
# subsetData then retains its database signaling genes, complexes and cofactors.
root<-normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/')
.libPaths(c(file.path(root,'runtime/R-library'),.libPaths()))
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(data.table);library(jsonlite)})
future::plan('sequential');options(future.globals.maxSize=4*1024^3)
cfg<-fromJSON(file.path(root,'cellchat_config.json'));set.seed(cfg$seed)
inp<-file.path(root,'data/CellChat');out<-file.path(root,'derived/CellChat');dir.create(out,showWarnings=FALSE)
jobs<-fread(file.path(inp,'jobs.csv'))
args<-commandArgs(trailingOnly=TRUE);if(length(args))jobs<-jobs[grepl(args[1],job)]
db<-CellChatDB.human
# Focus on protein-mediated intercellular signaling relevant to tumor immunity.
db$interaction<-db$interaction[db$interaction$annotation!='Non-protein Signaling',,drop=FALSE]
if(!file.exists(file.path(inp,'CellChatDB_protein_locked.rds')))saveRDS(db,file.path(inp,'CellChatDB_protein_locked.rds'))
if(!file.exists(file.path(root,'tables/cellchat_database_interactions.csv')))fwrite(as.data.table(db$interaction),file.path(root,'tables/cellchat_database_interactions.csv'))
for(k in seq_len(nrow(jobs))){
 j<-jobs[k];prefix<-file.path(inp,j$job);result<-file.path(out,paste0(j$job,'.rds'))
 if(file.exists(result)){cat('CHECKPOINT',j$job,'\n');next}
 stopifnot(file.exists(paste0(prefix,'_extraction.json')))
 cat('START',j$job,format(Sys.time()),'\n');flush.console()
 info<-fromJSON(paste0(prefix,'_extraction.json'));m<-fread(paste0(prefix,'_cells.csv'));setorder(m,local_col)
 genes<-readLines(paste0(prefix,'_genes.txt'));stopifnot(!anyDuplicated(genes))
 con<-file(paste0(prefix,'_counts.bin'),'rb');trip<-readBin(con,integer(),n=info$nnz*3,size=4,endian='little');close(con)
 stopifnot(length(trip)==info$nnz*3,all(trip>0));dim(trip)<-c(3,info$nnz)
 x<-sparseMatrix(i=trip[1,],j=trip[2,],x=trip[3,],dims=c(length(genes),nrow(m)),dimnames=list(genes,m$barcode));rm(trip)
 lib<-Matrix::colSums(x);stopifnot(all(lib>0))
 # Independently check the selected cells against the prior full-library extraction.
 refname<-if(j$dataset=='GEO')'GEO_GSE236581' else j$dataset
 ref<-fread(file.path(root,'derived',paste0(refname,'_expression.csv.gz')),select=c('barcode','library','METTL3','STRAP','PTBP1'))
 ix<-match(m$barcode,ref$barcode);stopifnot(!anyNA(ix),all(lib==ref$library[ix]))
 for(g in c('METTL3','STRAP','PTBP1'))stopifnot(all(as.numeric(x[g,])==ref[[g]][ix]))
 fwrite(data.table(job=j$job,cells=nrow(m),genes=nrow(x),all_libraries_and_targets_match=TRUE),file.path(root,'audit',paste0(j$job,'_cellchat_input_validation.csv')))
 x<-x%*%Diagonal(x=10000/lib);x@x<-log1p(x@x);colnames(x)<-m$barcode
 meta<-data.frame(labels=factor(m$label),samples=factor(j$sample),row.names=m$barcode)
 cc<-createCellChat(x,meta=meta,group.by='labels');cc@DB<-db;cc<-subsetData(cc)
 cc<-identifyOverExpressedGenes(cc,do.fast=FALSE)
 cc<-identifyOverExpressedInteractions(cc)
 cat('CANDIDATE_LR',nrow(cc@LR$LRsig),'signaling genes',nrow(cc@data.signaling),'\n');flush.console()
 cc<-computeCommunProb(cc,type=cfg$type,raw.use=TRUE,population.size=FALSE,nboot=cfg$nboot,seed.use=cfg$seed)
 cc<-filterCommunication(cc,min.cells=cfg$min_group_cells)
 cc<-computeCommunProbPathway(cc);cc<-aggregateNet(cc,remove.isolate=FALSE)
 stopifnot(all(is.finite(cc@net$prob)),all(cc@net$prob>=0),all(cc@net$pval>=0&cc@net$pval<=1))
 # Compact exact CellChat results retain all LR probabilities/p-values, pathway
 # arrays, labels, parameters, database identifiers and input validation.
 saveRDS(list(job=as.list(j),net=cc@net,netP=cc@netP,LR=cc@LR,options=cc@options,
              cell_counts=table(cc@idents),db_genes_present=rownames(cc@data.signaling),
              full_gene_library=lib,config=cfg,package_version=as.character(packageVersion('CellChat'))),result)
 cat('COMPLETE',j$job,'pathways',length(cc@netP$pathways),format(Sys.time()),'\n');flush.console()
 rm(cc,x,ref,meta);gc()
}
writeLines(capture.output(sessionInfo()),file.path(root,'audit/CellChat_final_sessionInfo.txt'))
