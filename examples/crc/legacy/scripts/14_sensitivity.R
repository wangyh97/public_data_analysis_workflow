# Sensitivity analyses on independent original studies, keeping all null results.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(ggplot2)})
root<-normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/')
cfg<-fromJSON(file.path(root,'analysis_config.json'));set.seed(cfg$seed)
d<-fread(file.path(root,'tables/tumor_patient_analysis.csv'))
a<-fread(file.path(root,'tables/sample_compartment_summary.csv'))
primary<-function(x)x[compartment=='Malignant'&time=='Pre'&((portal=='GEO'&study=='GSE236581')|(portal=='scCT-DB'&study=='GSE205506'))]
z<-primary(d);outcomes<-c('score_MHC_I_core','score_APM');loo<-list();sensitivity<-list();adjust<-list();highlow<-list();highlow_tests<-list()
for(st in unique(z$study))for(g in cfg$targets){
 v<-z[study==st];x<-paste0('mean_',g)
 med<-median(v[[x]]);v[,expression_group:=ifelse(get(x)>med,'High','Low')]
 highlow[[length(highlow)+1]]<-v[,.(portal,study,patient,gene=g,target_mean=get(x),median_cutoff=med,expression_group,n_cells,final_response,final_regimen)]
 for(y in outcomes){
  vv<-v[is.finite(get(x))&is.finite(get(y))]
  if(nrow(vv)>=6){
   vals<-vapply(seq_len(nrow(vv)),function(i)cor(vv[-i][[x]],vv[-i][[y]],method='spearman'),numeric(1));rho<-cor(vv[[x]],vv[[y]],method='spearman')
   loo[[length(loo)+1]]<-data.table(study=st,gene=g,outcome=y,n=nrow(vv),full_rho=rho,LOO_min=min(vals),LOO_max=max(vals),LOO_same_sign=all(sign(vals)==sign(rho)))
   vv[,combo:=as.integer(grepl('Capecitabine|Oxaliplatin|Fluorouracil|FOLFOX',final_regimen,ignore.case=TRUE))]
   for(group in c('ICI_only','Combination')){
    ss<-vv[combo==as.integer(group=='Combination')&!is.na(final_regimen)&final_regimen!='']
    if(nrow(ss)>=6&&sd(ss[[x]])>0&&sd(ss[[y]])>0){ct<-cor.test(ss[[x]],ss[[y]],method='spearman',exact=FALSE);sensitivity[[length(sensitivity)+1]]<-data.table(study=st,gene=g,outcome=y,analysis=group,n=nrow(ss),rho=unname(ct$estimate),p=ct$p.value)}
   }
   if(nrow(vv)>=10){
    cc<-data.frame(IFN=rank(vv$score_IFN),Cycle=rank(vv$score_Cycle),Complexity=rank(vv$mean_nFeature))
    if(uniqueN(vv$combo)>1)cc$Combination<-vv$combo
    rx<-resid(lm(rank(vv[[x]])~.,data=cc));ry<-resid(lm(rank(vv[[y]])~.,data=cc));r<-cor(rx,ry);df<-nrow(vv)-ncol(cc)-2
    pp<-2*pt(-abs(r*sqrt(df/(1-r*r))),df)
    adjust[[length(adjust)+1]]<-data.table(study=st,gene=g,outcome=y,n=nrow(vv),rho=r,p=pp,covariates=paste(names(cc),collapse=' + '))
   }
  }
 }
 for(y in intersect(c('T_cytotoxic','T_dysfunction','T_stem_memory','T_fraction_all','T_fraction_immune','CD8_cytotoxic','CD8_dysfunction','CD8_fraction_all'),names(v))){
  vv<-v[is.finite(get(y))];nn<-table(vv$expression_group)
  if(length(nn)==2&&min(nn)>=3){ct<-wilcox.test(vv[[y]]~vv$expression_group,exact=FALSE);highlow_tests[[length(highlow_tests)+1]]<-data.table(study=st,gene=g,outcome=y,n_low=nn['Low'],n_high=nn['High'],median_High_minus_Low=median(vv[expression_group=='High'][[y]])-median(vv[expression_group=='Low'][[y]]),p=ct$p.value)}
 }
}
# Change the tumor-cell eligibility threshold without redefining malignant cells.
for(cut in c(10,30,50,100)){
 aa<-primary(a)[tumor_tissue==TRUE&n_cells>=cut]
 gg<-c(paste0('mean_',cfg$targets),outcomes)
 aa<-aa[,lapply(.SD,mean),by=.(portal,study,patient),.SDcols=gg]
 for(st in unique(aa$study))for(g in cfg$targets)for(y in outcomes){
  vv<-aa[study==st];x<-paste0('mean_',g)
  if(nrow(vv)>=6&&sd(vv[[x]])>0&&sd(vv[[y]])>0){ct<-cor.test(vv[[x]],vv[[y]],method='spearman',exact=FALSE);sensitivity[[length(sensitivity)+1]]<-data.table(study=st,gene=g,outcome=y,analysis=paste0('min_tumor_cells_',cut),n=nrow(vv),rho=unname(ct$estimate),p=ct$p.value)}
 }
}
save_table<-function(x,name){r<-rbindlist(x,fill=TRUE);if(nrow(r)&&'p'%in%names(r))r[,q:=p.adjust(p,'BH'),by=study];fwrite(r,file.path(root,'tables',name));r}
lo<-save_table(loo,'leave_one_patient_out.csv');ss<-save_table(sensitivity,'regimen_and_cell_threshold_sensitivity.csv');ad<-save_table(adjust,'primary_IFN_cycle_complexity_regimen_adjusted.csv');hl<-save_table(highlow,'patient_high_low_groups.csv');ht<-save_table(highlow_tests,'patient_high_low_T_state_tests.csv')
fwrite(z,file.path(root,'tables/primary_independent_studies.csv'))
# Clear, baseline-only figures comparing original studies; portals are not replicates.
r<-fread(file.path(root,'tables/correlations_all.csv'));r<-r[time=='Pre'&compartment=='Malignant'&outcome%in%outcomes&((portal=='GEO'&study=='GSE236581')|(portal=='scCT-DB'&study=='GSE205506'))]
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(strip.background=element_blank(),legend.position='bottom'))
r[,gene:=factor(gene,levels=cfg$targets)]
p<-ggplot(r,aes(rho,gene,color=study))+geom_vline(xintercept=0,color='grey60',linetype=2)+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',position=position_dodge(width=.5),width=.15)+geom_point(position=position_dodge(width=.5),size=2)+facet_wrap(~outcome,ncol=1)+scale_color_manual(values=c(GSE236581='#287D8E',GSE205506='#CE7952'))+labs(x='Spearman rho (patient bootstrap 95% interval)',y=NULL,title='Baseline malignant-cell antigen-presentation associations',subtitle='Independent studies; full GEO GSE236581 and paired scCT-DB GSE205506')
ggsave(file.path(root,'figures/07_primary_MHC_independent_studies.pdf'),p,width=183,height=125,units='mm',device=cairo_pdf);ggsave(file.path(root,'figures/07_primary_MHC_independent_studies.png'),p,width=183,height=125,units='mm',dpi=300)
# Scatter panels show all patients, not only significant candidate/outcome pairs.
sc<-melt(z,id.vars=c('study','patient','final_response','score_MHC_I_core','score_APM'),measure.vars=paste0('mean_',cfg$targets),variable.name='gene',value.name='target')
sc[,gene:=sub('mean_','',gene)];sc<-melt(sc,id.vars=c('study','patient','final_response','gene','target'),measure.vars=outcomes,variable.name='signature',value.name='score')
fwrite(sc,file.path(root,'tables/primary_scatter_source.csv'))
p<-ggplot(sc,aes(target,score,color=final_response))+geom_point(size=1.7)+facet_grid(study+signature~gene,scales='free')+scale_color_manual(values=c(R='#287D8E',NR='#CE7952'))+labs(x='Target mean log1p(CP10K)',y='Mean signature expression',color='Response',title='Patient-level source data for primary associations')
ggsave(file.path(root,'figures/08_primary_scatter.pdf'),p,width=183,height=180,units='mm',device=cairo_pdf);ggsave(file.path(root,'figures/08_primary_scatter.png'),p,width=183,height=180,units='mm',dpi=300)
writeLines(capture.output(sessionInfo()),file.path(root,'audit/R_final_sessionInfo.txt'))
cat('SENSITIVITY COMPLETE',nrow(z),'primary patient rows\n')
