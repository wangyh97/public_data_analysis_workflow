# Two-sided PAIRED t tests on patient pseudobulk log2(CPM+1).
# Each comparison matches patient IDs explicitly. No unpaired or cell-level tests.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
this_script<-sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])
tt_panel<-normalizePath(file.path(dirname(this_script),'..'),winslash='/')
tt_d<-fread(file.path(tt_panel,'source/trajectory_plot_data.csv'))
tt_pre<-tt_d[stage=='Pre',.(cohort,patient,gene,pre=expression)]
tt_pairs<-merge(tt_d[stage!='Pre'],tt_pre,by=c('cohort','patient','gene'),all=FALSE)
stopifnot(!anyDuplicated(tt_pairs[,.(cohort,patient,gene,stage)]))
tt_pairs[,difference:=expression-pre]
tt_result<-rbindlist(lapply(c('All','R','NR'),function(gr){
 z<-if(gr=='All')tt_pairs else tt_pairs[response==gr]
 z[,{
  n<-.N;se<-if(n>=2)sd(difference)/sqrt(n) else NA_real_
  valid<-n>=2 && is.finite(se) && se>0
  test<-if(valid)t.test(expression,pre,paired=TRUE,alternative='two.sided',conf.level=.95) else NULL
  list(group=gr,n_pairs=n,n_R=sum(response=='R'),n_NR=sum(response=='NR'),mean_pre=mean(pre),mean_post=mean(expression),mean_difference=mean(difference),
   sd_difference=if(n>=2)sd(difference) else NA_real_,t=if(valid)unname(test$statistic) else NA_real_,df=n-1,
   ci_low=if(valid)test$conf.int[1] else NA_real_,ci_high=if(valid)test$conf.int[2] else NA_real_,
   p=if(valid)test$p.value else NA_real_,status=if(n<2)'not testable: fewer than 2 pairs' else if(!valid)'not testable: zero variance of differences' else if(n==2)'tested; n=2, highly unstable' else 'tested')
 },by=.(cohort,gene,stage)]
}),fill=TRUE)
tt_result[,q:=p.adjust(p,'BH'),by=.(cohort,group)]
fwrite(tt_pairs,file.path(tt_panel,'source/paired_t_test_pairs.csv'))
fwrite(tt_result,file.path(tt_panel,'source/paired_t_test_statistics.csv'))
print(tt_result[group=='All',.(cohort,gene,stage,n_pairs,mean_difference,p,q,status)])
