# Run from any working directory: Rscript data_processing.R [--rebuild-raw]
# Same-tumor-sample malignant target expression is paired with All_T or CD8_T
# module expression BEFORE equal-sample-weight patient aggregation. All_T and
# CD8_T overlap and are not independent replicates. See config.json for genes.
script<-sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])
panel_dir<-normalizePath(file.path(dirname(script),'..'),winslash='/')
source(file.path(panel_dir,'../_shared/process.R'))
process_panel(panel_dir,'--rebuild-raw'%in%commandArgs(TRUE))
