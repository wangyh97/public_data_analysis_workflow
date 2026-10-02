# Match original GEO cell annotations with the public clinical sample tables.
# Do not guess response or CRC diagnosis for patients missing from the CRC tables.
suppressPackageStartupMessages(library(data.table))
root<-normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/')
m<-as.data.table(read.table(gzfile(file.path(root,'data/GEO/GSE236581_CRC-ICB_metadata.txt.gz')),header=TRUE,row.names=1),keep.rownames='barcode')
s<-fread(file.path(root,'audit/official_sample_tables.csv'))
c<-unique(s[grepl('GSE236581',dataset),.(`Patient ID`,Response)])
stopifnot(!anyDuplicated(c[['Patient ID']]))
ref<-fread(file.path(root,'derived/GSE236581_Colorectal cancer_1_metadata.csv.gz'),select=c('patient','drug_name'))
reg<-unique(ref);stopifnot(!anyDuplicated(reg$patient))
m[,`:=`(patient=Patient,sample_raw=Ident,Sample_Location=Tissue,Response=c$Response[match(Patient,c[['Patient ID']])],drug_name=reg$drug_name[match(Patient,reg$patient)],is_crc=Patient%in%c[['Patient ID']])]
m[,cluster:=fcase(MajorCellType=='B','B Cell',MajorCellType=='Epi','Epithelia Cell',MajorCellType=='ILC','Innate Lymphocyte',MajorCellType=='Mye','Myeloid Cell',MajorCellType=='Stromal','Stromal Cell',MajorCellType=='T','T Cell',default='Unknown')]
m[,treatment:=ifelse(Treatment=='I','Pre','Post')]
fwrite(m,file.path(root,'derived/GEO_GSE236581_metadata.csv.gz'))
fwrite(m[,.(n_cells=.N,n_samples=uniqueN(Ident)),by=.(Patient,is_crc,Response,drug_name)],file.path(root,'audit/GEO236581_patient_clinical_mapping.csv'))
cat('Author metadata cells',nrow(m),'CRC-mapped cells',sum(m$is_crc),'CRC patients',uniqueN(m[is_crc==TRUE]$patient),'\n')
