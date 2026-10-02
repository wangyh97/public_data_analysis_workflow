# Run from any working directory: Rscript data_processing.R [--rebuild-raw]
# Counts -> full-library log1pCP10K -> sample means -> paired patient means
# -> Spearman correlations. Configuration: config.json; raw input paths are
# resolved relative to CRC_single_cell. The shared engine is distributed with
# this panel bundle to keep the two panels' normalization identical.
script<-sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])
panel_dir<-normalizePath(file.path(dirname(script),'..'),winslash='/')
source(file.path(panel_dir,'../_shared/process.R'))
process_panel(panel_dir,'--rebuild-raw'%in%commandArgs(TRUE))
