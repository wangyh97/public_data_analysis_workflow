# User requested completion of the missing CellChat analysis. Install only into
# this project's dedicated library; existing user/system libraries are unchanged.
root <- normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/')
lib <- file.path(root,'runtime/R-library');dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()));options(timeout=600,repos=c(CRAN='https://cloud.r-project.org'))
dir.create(file.path(root,'runtime/downloads'),recursive=TRUE,showWarnings=FALSE)
manifest <- jsonlite::fromJSON('https://blaserlab.r-universe.dev/api/packages/CellChat')
jsonlite::write_json(manifest,file.path(root,'audit/cellchat_binary_provenance.json'),pretty=TRUE,auto_unbox=TRUE)
stopifnot(manifest$RemoteSha=='75253cd0c9e68410e6e721a6d3a0419a1d7e358f',manifest$RemoteUrl=='https://github.com/jinworks/CellChat')
deps <- manifest$`_dependencies`;needed<-unique(deps$package[deps$role %in% c('Depends','Imports','LinkingTo')]);needed<-setdiff(needed,'R')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
writeLines(c(paste('Dedicated library:',lib),paste('Missing required packages:',paste(missing,collapse=', ')),
 'CRAN dependencies: Windows R 4.4 binaries; BiocNeighbors: Bioconductor 3.20 binary.',
 'CellChat: r-universe build of official jinworks source, commit 75253cd; R 4.5 build requires explicit runtime smoke test on R 4.4.'),file.path(root,'audit/cellchat_install_plan.txt'))
cran <- setdiff(missing,c('BiocNeighbors','ComplexHeatmap','BiocGenerics'))
if(length(cran))install.packages(cran,lib=lib,type='win.binary',destdir=file.path(root,'runtime/downloads'))
if(!requireNamespace('BiocNeighbors',quietly=TRUE))BiocManager::install('BiocNeighbors',lib=lib,ask=FALSE,update=FALSE,type='binary')
if(!requireNamespace('CellChat',quietly=TRUE)){
 b<-manifest$`_binaries`;u<-b$fileid[b$os=='win' & b$arch=='x86_64' & b$r=='4.5.3'];stopifnot(length(u)==1)
 dest<-file.path(root,'runtime/downloads/CellChat_2.2.0.9001.zip');download.file(u,dest,mode='wb')
 install.packages(dest,repos=NULL,type='win.binary',lib=lib)
}
library(CellChat)
stopifnot(as.character(packageVersion('CellChat'))=='2.2.0.9001')
writeLines(capture.output(sessionInfo()),file.path(root,'audit/CellChat_sessionInfo.txt'))
cat('CELLCHAT LOAD PASSED\n')
