# Patient-level inference and figures. All eligible tests are retained, including
# nonsignificant results. No cell-level p-values are used as patient evidence.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite)})
root<-normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/');cfg<-fromJSON(file.path(root,'analysis_config.json'));set.seed(cfg$seed)
a<-fread(file.path(root,'tables/sample_compartment_summary.csv'));targets<-cfg$targets
save_plot<-function(p,name,w=183,h=125){
  ggsave(file.path(root,'figures',paste0(name,'.pdf')),p,width=w,height=h,units='mm',device=cairo_pdf)
  ggsave(file.path(root,'figures',paste0(name,'.png')),p,width=w,height=h,units='mm',dpi=300)
}
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=8),legend.position='bottom',plot.title=element_text(size=11,face='bold'),plot.subtitle=element_text(size=8)))
colors<-c(R='#287D8E',NR='#CE7952',Pre='#287D8E',Post='#CE7952')
# Composition is a proportion of recovered cells, not an absolute infiltration density.
comp<-a[compartment!='CD8_T',.(n_all=sum(n_cells),n_T=sum(n_cells[compartment=='T']),n_immune=sum(n_cells[compartment%in%c('T','NK','NK_ILC','B_plasma','Myeloid')])),by=.(portal,study,sample,patient,time)]
comp[,`:=`(T_fraction_all=n_T/pmax(n_all,1),T_fraction_immune=n_T/pmax(n_immune,1))]
tstate<-a[compartment=='T'&n_cells>=cfg$min_T_cells,.(portal,study,sample,patient,time,T_cytotoxic=score_Cytotoxic,T_dysfunction=score_Dysfunction,T_stem_memory=score_Stem_memory,T_Treg_program=score_Treg)]
cd8state<-a[compartment=='CD8_T'&n_cells>=cfg$min_T_cells,.(portal,study,sample,patient,time,CD8_cytotoxic=score_Cytotoxic,CD8_dysfunction=score_Dysfunction,CD8_stem_memory=score_Stem_memory,n_CD8=n_cells)]
cd8counts<-a[compartment=='CD8_T',.(portal,study,sample,patient,time,n_CD8_total=n_cells)]
cd8_annotated_studies<-unique(paste(cd8counts$portal,cd8counts$study))
tum<-a[tumor_tissue==TRUE&compartment%in%c('Malignant','Epithelial_uncertain')&n_cells>=cfg$min_tumor_cells]
tum<-merge(tum,comp,by=c('portal','study','sample','patient','time'),all.x=TRUE)
tum<-merge(tum,tstate,by=c('portal','study','sample','patient','time'),all.x=TRUE)
tum<-merge(tum,cd8state,by=c('portal','study','sample','patient','time'),all.x=TRUE)
tum<-merge(tum,cd8counts,by=c('portal','study','sample','patient','time'),all.x=TRUE)
tum[is.na(n_CD8_total)&paste(portal,study)%in%cd8_annotated_studies,n_CD8_total:=0]
tum[,CD8_fraction_all:=n_CD8_total/n_all]
tum[n_all<cfg$min_sample_cells_for_composition,c('T_fraction_all','T_fraction_immune','CD8_fraction_all'):=.(NA_real_,NA_real_,NA_real_)]
key<-c('portal','study','patient','time','compartment','final_response','final_regimen')
num<-names(tum)[vapply(tum,is.numeric,logical(1))]
# Multiple tumor specimens at one timepoint receive equal sample weight here.
# Samples remain separately available for inspection; patients contribute once.
d<-tum[,lapply(.SD,function(v)if(all(is.na(v)))NA_real_ else mean(v,na.rm=TRUE)),by=key,.SDcols=num]
fwrite(tum,file.path(root,'tables/tumor_sample_analysis.csv'));fwrite(d,file.path(root,'tables/tumor_patient_analysis.csv'));fwrite(comp,file.path(root,'tables/cell_composition.csv'))
outcomes<-c('score_MHC_I_core','score_APM','mean_HLA-A','mean_HLA-B','mean_HLA-C','mean_B2M','T_cytotoxic','T_dysfunction','T_stem_memory','T_Treg_program','T_fraction_all','T_fraction_immune','CD8_cytotoxic','CD8_dysfunction','CD8_stem_memory','CD8_fraction_all')
strata<-unique(d[,.(portal,study,time,compartment)])
results<-list();adjusted<-list();response<-list();eligibility<-list()
bootrho<-function(x,y){
  v<-replicate(cfg$bootstrap_replicates,{i<-sample.int(length(x),replace=TRUE);if(sd(x[i])==0||sd(y[i])==0)NA_real_ else cor(x[i],y[i],method='spearman')})
  if(sum(is.finite(v))<cfg$bootstrap_replicates*.8)c(NA_real_,NA_real_) else quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE)
}
for(k in seq_len(nrow(strata))){
 s<-strata[k];z<-d[portal==s$portal&study==s$study&time==s$time&compartment==s$compartment]
 nr<-sum(z$final_response=='NR',na.rm=TRUE);rr<-sum(z$final_response=='R',na.rm=TRUE)
 eligibility[[k]]<-cbind(s,data.table(n_patients=nrow(z),n_R=rr,n_NR=nr,response_test_eligible=min(nr,rr)>=cfg$min_patients_per_response_group))
 for(g in targets){
  xname<-paste0('mean_',g)
  if(min(nr,rr)>=cfg$min_patients_per_response_group){
   v<-z[final_response%in%c('R','NR')];test<-wilcox.test(v[[paste0('pb_',g)]]~v$final_response,exact=FALSE)
   response[[length(response)+1]]<-cbind(s,data.table(gene=g,n_R=rr,n_NR=nr,median_R_minus_NR=median(v[final_response=='R'][[paste0('pb_',g)]])-median(v[final_response=='NR'][[paste0('pb_',g)]]),p=test$p.value))
  }
  for(yname in outcomes){
   v<-z[is.finite(get(xname))&is.finite(get(yname))];n<-nrow(v)
   if(n<cfg$min_patients_correlation||sd(v[[xname]])==0||sd(v[[yname]])==0)next
   ct<-cor.test(v[[xname]],v[[yname]],method='spearman',exact=FALSE);ci<-bootrho(v[[xname]],v[[yname]])
   results[[length(results)+1]]<-cbind(s,data.table(gene=g,outcome=yname,n=n,rho=unname(ct$estimate),ci_low=ci[1],ci_high=ci[2],p=ct$p.value))
   if(yname%in%c('score_MHC_I_core','score_APM')&&n>=cfg$min_patients_adjusted_correlation){
    # Partial rank correlation: remove shared IFN/cycle/complexity variation.
    cv<-c('score_IFN','score_Cycle','mean_nFeature');v<-v[complete.cases(v[,..cv])]
    if(nrow(v)>=cfg$min_patients_adjusted_correlation){
     cc<-as.data.frame(lapply(v[,..cv],rank));xx<-rank(v[[xname]]);yy<-rank(v[[yname]])
     rx<-resid(lm(xx~.,data=cc));ry<-resid(lm(yy~.,data=cc));rho<-cor(rx,ry);df<-nrow(v)-ncol(cc)-2
     pp<-2*pt(-abs(rho*sqrt(df/(1-rho^2))),df=df)
     adjusted[[length(adjusted)+1]]<-cbind(s,data.table(gene=g,outcome=yname,n=nrow(v),rho=rho,p=pp,covariates=paste(cv,collapse=' + ')))
    }
   }
  }
 }
}
res<-rbindlist(results,fill=TRUE);adj<-rbindlist(adjusted,fill=TRUE);rsp<-rbindlist(response,fill=TRUE)
if(nrow(res)){
 res[,family:=ifelse(time=='Pre'&compartment=='Malignant'&outcome%in%c('score_MHC_I_core','score_APM'),'primary_MHC','exploratory')]
 res[,q:=p.adjust(p,'BH'),by=.(portal,study,family)]
}
if(nrow(adj))adj[,q:=p.adjust(p,'BH'),by=.(portal,study)]
if(nrow(rsp))rsp[,q:=p.adjust(p,'BH'),by=.(portal,study,time,compartment)]
fwrite(res,file.path(root,'tables/correlations_all.csv'));fwrite(adj,file.path(root,'tables/correlations_IFN_adjusted.csv'));fwrite(rsp,file.path(root,'tables/response_tests.csv'));fwrite(rbindlist(eligibility),file.path(root,'tables/analysis_eligibility.csv'))
# Descriptive all-cell-type target expression: average across patient/time samples.
dot<-melt(a,id.vars=c('portal','study','patient','time','compartment'),measure.vars=paste0('mean_',targets),variable.name='gene',value.name='expression')
dot[,gene:=sub('mean_','',gene)];dot<-dot[,.(expression=mean(expression,na.rm=TRUE)),by=.(portal,study,patient,time,compartment,gene)]
dot<-dot[,.(expression=mean(expression),n_patients=uniqueN(patient)),by=.(portal,study,compartment,gene)]
fwrite(dot,file.path(root,'tables/figure_expression_source.csv'))
p<-ggplot(dot,aes(gene,compartment,color=expression,size=n_patients))+geom_point()+facet_wrap(~portal+study,ncol=3)+scale_color_viridis_c(option='C')+scale_size(range=c(1.5,5))+labs(x=NULL,y=NULL,color='Mean log1p(CP10K)',size='Patients',title='Target expression across cell compartments')
save_plot(p,'01_target_expression',183,140)
if(nrow(d)){
 v<-melt(d,id.vars=key,measure.vars=paste0('pb_',targets),variable.name='gene',value.name='expression');v[,gene:=sub('pb_','',gene)]
 fwrite(v,file.path(root,'tables/figure_response_source.csv'))
 for(port in unique(v$portal))for(st in unique(v[portal==port]$study)){
  available<-v[portal==port&study==st&final_response%in%c('R','NR')]
  for(tm in unique(available$time)){
  vv<-available[time==tm]
  p<-ggplot(vv,aes(final_response,expression,color=final_response))+geom_boxplot(width=.5,outlier.shape=NA,alpha=.2)+geom_point(position=position_jitter(width=.10,seed=cfg$seed),size=1.8)+facet_grid(time+compartment~gene,scales='free_y')+scale_color_manual(values=colors)+labs(x='Clinical response',y='Target log2(pseudobulk CPM + 1)',title=paste(st,port),subtitle='One point per patient/time; group sample sizes and test eligibility are in analysis_eligibility.csv')
  save_plot(p,paste0('02_response_',port,'_',st,if(tm=='Pre')'' else paste0('_',tm)),183,105)
  }
 }
}
if(nrow(res)){
 for(port in unique(res$portal))for(st in unique(res[portal==port]$study)){
  for(tm in unique(res[portal==port&study==st]$time)){
  r<-res[portal==port&study==st&time==tm];r[,gene:=factor(gene,levels=targets)];r[,label:=sprintf('%.2f%s\nn=%d',rho,ifelse(q<.05,'*',''),n)]
  p<-ggplot(r,aes(gene,outcome,fill=rho))+geom_tile(color='white')+geom_text(aes(label=label),size=2.3)+facet_wrap(~time+compartment,nrow=1)+scale_fill_gradient2(low='#3B6FA0',mid='white',high='#BA5C4F',limits=c(-1,1))+labs(x=NULL,y=NULL,title=paste('Patient-level correlations:',st,port),subtitle='Spearman rho; * BH q < 0.05. n = independent patients.',fill='rho')
  save_plot(p,paste0('03_correlations_',port,'_',st,if(tm=='Pre')'' else paste0('_',tm)),183,155)
  }
 }
 # Show every qualifying MHC test, avoiding selection by statistical significance.
 r<-res[time=='Pre'&outcome%in%c('score_MHC_I_core','score_APM')];r[,label:=paste(portal,study,time,compartment,sep=' | ')]
 p<-ggplot(r,aes(rho,label,color=gene))+geom_vline(xintercept=0,linetype=2,color='grey60')+geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation='y',position=position_dodge(width=.6),width=.2)+geom_point(position=position_dodge(width=.6),size=1.6)+facet_wrap(~outcome,ncol=1)+scale_color_manual(values=c(METTL3='#287D8E',STRAP='#CE7952',PTBP1='#75618B'))+labs(x='Spearman rho (patient bootstrap 95% interval)',y=NULL,title='MHC-I and antigen-processing associations')
 save_plot(p,'04_MHC_forest',183,max(130,uniqueN(r$label)*9))
}
# Within-sample cell correlations summarized per patient, with no pseudoreplication.
w<-fread(file.path(root,'tables/within_sample_correlations.csv'))
if(nrow(w)){
 w<-w[tumor_tissue==TRUE,.(rho=mean(rho)),by=.(portal,study,patient,time,compartment,gene,signature)]
 fwrite(w,file.path(root,'tables/within_patient_correlations.csv'))
 p<-ggplot(w[time=='Pre'],aes(gene,rho,color=gene))+geom_hline(yintercept=0,color='grey60',linetype=2)+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=.12,seed=cfg$seed),size=.8,alpha=.6)+facet_grid(portal+study+compartment~signature)+labs(x=NULL,y='Within-sample Spearman rho',title='Cell-level association consistency across patients',subtitle='Baseline only; each point is a patient/compartment. Descriptive, no cell-level p-values.')+guides(color='none')
 save_plot(p,'05_within_patient_MHC',183,210)
}
# Paired baseline/post target change where the same patient has enough tumor cells.
paired<-list()
for(g in targets){
 z<-d[time%in%c('Pre','Post','Post_II'),c('portal','study','patient','time','compartment','final_response',paste0('pb_',g)),with=FALSE];z[time=='Post_II',time:='Post'];setnames(z,paste0('pb_',g),'expression')
 wide<-dcast(z,portal+study+patient+compartment+final_response~time,value.var='expression')
 if(all(c('Pre','Post')%in%names(wide))){wide<-wide[is.finite(Pre)&is.finite(Post)];wide[,`:=`(gene=g,delta_Post_minus_Pre=Post-Pre)];paired[[g]]<-wide}
}
pa<-rbindlist(paired,fill=TRUE);fwrite(pa,file.path(root,'tables/paired_changes.csv'))
if(nrow(pa)){
 p<-ggplot(pa,aes(gene,delta_Post_minus_Pre,color=final_response))+geom_hline(yintercept=0,linetype=2,color='grey60')+geom_point(position=position_jitter(width=.12,seed=cfg$seed),size=1.5)+facet_wrap(~portal+study+compartment)+scale_color_manual(values=colors,na.value='grey60')+labs(x=NULL,y='Paired change in log2(pseudobulk CPM + 1)',title='Treatment-associated target changes',subtitle='One point per paired patient; response association is not treatment causality')
 save_plot(p,'06_paired_changes',183,125)
}
cat('FIGURES COMPLETE; patient rows',nrow(d),'correlation tests',nrow(res),'response tests',nrow(rsp),'\n')
