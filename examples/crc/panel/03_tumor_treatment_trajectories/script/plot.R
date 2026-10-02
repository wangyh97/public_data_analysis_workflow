# Fig3: R-only patient pseudobulk trajectories. Every point is an observed node.
# Plotting is independent of raw matrices; source/trajectory_plot_data.csv is enough.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite)})
script<-sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])
panel<-normalizePath(file.path(dirname(script),'..'),winslash='/')
cfg<-fromJSON(file.path(panel,'script/config.json'))
# ================= STYLE / 用户自定义接口 ============================
STYLE<-list(
 font='Arial',font_size=9,title_size=11,strip_size=9,
 response_colors=c(R='#D73027',NR='#2166AC'), # responder red; non-responder blue
 line_width=.5,line_alpha=.65,point_size=1.9,point_alpha=.95,
 point_shapes=c(R=16,NR=17), # Shape reinforces the requested color distinction.
 overview_width_mm=183,overview_height_mm=135,
 patient_width_mm=183,patient_columns=4,patient_row_height_mm=39,
 patient_shared_y=TRUE, # FALSE permits separate y ranges per patient; not default.
 gene_order=c('METTL3','STRAP','PTBP1'),cohort_order=c('GSE236581','GSE205506'),
 dpi=600,formats=c('pdf','png'),show_overview=TRUE,show_patient_panels=TRUE,
 show_missing_bridges=TRUE, # dashed segment crosses an unobserved intermediate node
 x_padding=.13,
 show_t_tests=TRUE,stat_font_size=7.5,
 show_paired_boxes=TRUE,box_width_mm=183,box_height_mm=110,
 show_multistage_boxes=TRUE,multistage_height_mm=145,
 box_width=.48,jitter_width=.075,jitter_seed=20260910
)
# Change the values above and rerun. No pseudobulk/statistical recomputation needed.
# ================= END STYLE ========================================
d<-fread(file.path(panel,'source/trajectory_plot_data.csv'))
source(file.path(panel,'script/paired_t_tests.R'))
tt<-fread(file.path(panel,'source/paired_t_test_statistics.csv'))
paired_points<-fread(file.path(panel,'source/paired_t_test_pairs.csv'))
fmt_p<-function(v)ifelse(is.na(v),'NA',formatC(v,format='g',digits=3))
d[,gene:=factor(gene,levels=STYLE$gene_order)]
d[,response:=factor(response,levels=c('R','NR'))]
d[,patient_short:=sub('^GSE[0-9]+-','',patient)]
d[,patient_label:=paste0(patient_short,'  |  ',response)]
setorder(d,cohort,gene,patient,stage_order)
segments<-d[,.(x=stage_order,y=expression,xend=shift(stage_order,type='lead'),
 yend=shift(expression,type='lead'),response=response,patient_label=patient_label),by=.(cohort,gene,patient)]
segments<-segments[!is.na(xend)]
segments[,interval:=ifelse(xend-x>1,'Missing intermediate node','Adjacent observed nodes')]
if(!STYLE$show_missing_bridges)segments<-segments[xend-x==1]
theme_set(theme_classic(base_size=STYLE$font_size,base_family=STYLE$font)+theme(
 plot.title=element_text(size=STYLE$title_size,face='bold',margin=margin(b=6)),
 plot.subtitle=element_text(size=8,margin=margin(b=8)),
 plot.caption=element_text(size=7,hjust=0,margin=margin(t=8)),
 strip.background=element_blank(),strip.text=element_text(size=STYLE$strip_size,face='bold'),
 axis.text=element_text(color='#303030'),axis.line=element_line(linewidth=.35),
 axis.ticks=element_line(linewidth=.35),panel.spacing=grid::unit(6,'mm'),
 plot.margin=margin(9,9,9,9),legend.position='bottom',legend.title=element_blank()))
save_fig<-function(p,name,w,h){
 for(ext in STYLE$formats){
  dest<-file.path(panel,'figure',paste0(name,'.',ext))
  if(ext=='pdf'){grDevices::cairo_pdf(dest,width=w/25.4,height=h/25.4,family=STYLE$font);print(p);dev.off()}
  else if(ext=='png')ggsave(dest,p,width=w,height=h,units='mm',dpi=STYLE$dpi,device=grDevices::png,type='cairo')
  # Optional .svg support (requires svglite); default .pdf / .png use base R.
  else if(ext=='svg'){if(!requireNamespace('svglite',quietly=TRUE))stop('Install svglite to enable SVG');svglite::svglite(dest,width=w/25.4,height=h/25.4);print(p);dev.off()}
  else if(ext=='tiff')ggsave(dest,p,width=w,height=h,units='mm',dpi=STYLE$dpi,device=grDevices::tiff,type='cairo',compression='lzw')
  else stop('Supported formats: pdf, png, tiff, svg')
 }
}
draw<-function(dd,ss,stages){
 ggplot(dd,aes(stage_order,expression,color=response,shape=response))+
  geom_segment(data=ss,aes(x=x,y=y,xend=xend,yend=yend,color=response,linetype=interval),inherit.aes=FALSE,
    linewidth=STYLE$line_width,alpha=STYLE$line_alpha)+
  geom_point(size=STYLE$point_size,alpha=STYLE$point_alpha)+
  scale_color_manual(values=STYLE$response_colors,drop=FALSE,labels=c(R='Responder',NR='Non-responder'),limits=c('R','NR'))+
  scale_shape_manual(values=STYLE$point_shapes,drop=FALSE,labels=c(R='Responder',NR='Non-responder'),limits=c('R','NR'))+
  scale_linetype_manual(values=c('Adjacent observed nodes'='solid','Missing intermediate node'='dashed'),guide='none')+
  scale_x_continuous(breaks=seq_along(stages),labels=stages,limits=c(1, length(stages)),expand=expansion(add=STYLE$x_padding))+
  scale_y_continuous(expand=expansion(mult=c(.06,.10)))+
  labs(x='Treatment node',y='Tumor pseudobulk log2(CPM + 1)')
}
for(co in STYLE$cohort_order){
 dd<-d[cohort==co];ss<-segments[cohort==co];stages<-cfg$timepoints[[co]]
 pp<-unique(dd[,.(patient,response)])
 nlabel<-paste0(nrow(pp),' paired patients: ',sum(pp$response=='R'),' R / ',sum(pp$response=='NR'),' NR')
 endpoint<-if(co=='GSE236581')'RECIST response' else 'Pathological response (pCR / non-pCR)'
 note<-if(co=='GSE236581')'Author malignant annotation. Dashed lines bridge missing nodes; no values are imputed.' else 'scCT-DB malignant labels are unverified here; pCR Post cells require particular caution.'
 if(STYLE$show_overview){
  lab<-setNames(STYLE$gene_order,STYLE$gene_order)
  if(STYLE$show_t_tests)for(g in STYLE$gene_order){
   trow<-tt[cohort==co & gene==g & group=='All'];trow[,ord:=match(stage,stages)];setorder(trow,ord)
   lab[g]<-paste(c(g,paste0(trow$stage,' vs Pre | n=',trow$n_pairs,'\np=',fmt_p(trow$p),'; q=',fmt_p(trow$q))),collapse='\n')
  }
  p<-draw(dd,ss,stages)+facet_wrap(~gene,nrow=1,scales='free_y',labeller=labeller(gene=lab))+
   labs(title=paste0(co,' | Tumor expression during treatment'),subtitle=paste(nlabel,endpoint,sep='  |  '),
    caption=paste0('Two-sided paired t tests vs Pre; BH q across genes and follow-up nodes within cohort.\n',note))+
   theme(strip.text=element_text(size=STYLE$stat_font_size,lineheight=1.1))
  save_fig(p,paste0('Fig3_',co,'_trajectories'),STYLE$overview_width_mm,STYLE$overview_height_mm)
 }
 # 并列箱线图：每节点保留可用患者；检验仍按各随访节点与 Pre 单独配对。
 # x 轴 n 是箱体患者数，表头 paired n 是该项 t 检验的配对人数。
 if(STYLE$show_multistage_boxes && length(stages)>2){
  set.seed(STYLE$jitter_seed)
  offsets<-data.table(patient=sort(unique(dd$patient)))
  offsets[,offset:=runif(.N,-STYLE$jitter_width,STYLE$jitter_width)]
  bd<-merge(copy(dd),offsets,by='patient')
  bs<-merge(copy(ss),offsets,by='patient')
  counts<-unique(bd[,.(patient,stage_order)])[,.(n=.N),by=stage_order]
  labels<-paste0(stages,'\nn=',counts$n[match(seq_along(stages),counts$stage_order)])
  lab<-setNames(STYLE$gene_order,STYLE$gene_order)
  if(STYLE$show_t_tests)for(g in STYLE$gene_order){
   tx<-copy(tt[cohort==co & gene==g & group=='All'])
   tx[,ord:=match(stage,stages)];setorder(tx,ord)
   lab[g]<-paste(c(g,paste0(tx$stage,' vs Pre | paired n=',tx$n_pairs,
     '\np=',fmt_p(tx$p),'; q=',fmt_p(tx$q))),collapse='\n')
  }
  p<-ggplot(bd,aes(stage_order,expression))+
   geom_boxplot(aes(group=stage_order),width=STYLE$box_width,fill='#EFEFEF',color='#555555',linewidth=.4,outlier.shape=NA)+
   geom_segment(data=bs,aes(x=x+offset,xend=xend+offset,y=y,yend=yend,color=response,linetype=interval),inherit.aes=FALSE,linewidth=STYLE$line_width,alpha=STYLE$line_alpha)+
   geom_point(aes(x=stage_order+offset,color=response,shape=response),size=STYLE$point_size,alpha=STYLE$point_alpha,show.legend=TRUE)+
   facet_wrap(~gene,nrow=1,scales='free_y',labeller=labeller(gene=lab))+
   scale_x_continuous(breaks=seq_along(stages),labels=labels,limits=c(.6,length(stages)+.4))+
   scale_y_continuous(expand=expansion(mult=c(.07,.10)))+
   scale_color_manual(values=STYLE$response_colors,drop=FALSE,labels=c(R='Responder',NR='Non-responder'),limits=c('R','NR'))+
   scale_shape_manual(values=STYLE$point_shapes,drop=FALSE,labels=c(R='Responder',NR='Non-responder'),limits=c('R','NR'))+
   scale_linetype_manual(values=c('Adjacent observed nodes'='solid','Missing intermediate node'='dashed'),guide='none')+
   labs(title=paste0(co,' | Tumor expression across treatment nodes'),
    subtitle=paste0(nlabel,'; available observations at each node'),
    x='Treatment node',y='Tumor pseudobulk log2(CPM + 1)',
    caption=paste0('Box: median and IQR; whiskers: 1.5 IQR. Points and lines identify paired observations.\n',
      'Two-sided paired t vs Pre; BH q across 9 tests. Axis n: observed patients; paired n: matched patients.\n',
      'Dashed lines bridge missing nodes. Stage IV has only 2 patients; estimates are unstable.'))+
   theme(strip.text=element_text(size=STYLE$stat_font_size,lineheight=1.1))
  save_fig(p,paste0('Fig3_',co,'_all_stages_boxplot'),STYLE$box_width_mm,STYLE$multistage_height_mm)
 }
 if(STYLE$show_patient_panels)for(g in STYLE$gene_order){
  dg<-dd[gene==g];sg<-ss[gene==g]
  # Facets prevent same-color patient trajectories from becoming ambiguous.
  p<-draw(dg,sg,stages)+facet_wrap(~patient_label,ncol=STYLE$patient_columns,scales=if(STYLE$patient_shared_y)'fixed' else 'free_y',axes='all',axis.labels='all')+
   labs(title=paste0(co,' | ',g,' by patient'),subtitle=paste(nlabel,endpoint,sep='  |  '),
    caption=paste0('Each panel is one patient; shared y scale by default. Missing nodes have no points.\n',note))
  h<-ceiling(nrow(pp)/STYLE$patient_columns)*STYLE$patient_row_height_mm+45
  save_fig(p,paste0('Fig3_',co,'_',g,'_patients'),STYLE$patient_width_mm,h)
 }
 # Separate matched box+jitter comparison for each follow-up node. Every box
 # uses exactly the same patients on its two sides; Pre is reselected per node.
 if(STYLE$show_paired_boxes)for(tm in stages[-1]){
  pairs<-copy(paired_points[cohort==co & stage==tm])
  if(!nrow(pairs))next
  set.seed(STYLE$jitter_seed)
  offsets<-data.table(patient=sort(unique(pairs$patient)))
  offsets[,offset:=runif(.N,-STYLE$jitter_width,STYLE$jitter_width)]
  pairs<-merge(pairs,offsets,by='patient')
  boxes<-rbind(pairs[,.(patient,gene,response,node='Pre',position=1,expression=pre,offset)],
               pairs[,.(patient,gene,response,node=tm,position=2,expression=expression,offset)])
  boxes[,gene:=factor(gene,levels=STYLE$gene_order)]
  pairs[,gene:=factor(gene,levels=STYLE$gene_order)]
  tx<-tt[cohort==co & stage==tm & group=='All']
  tx[,gene:=factor(gene,levels=STYLE$gene_order)]
  ranges<-boxes[,.(ymin=min(expression),ymax=max(expression)),by=gene]
  tx<-merge(tx,ranges,by='gene');tx[,pad:=pmax(ymax-ymin,.2)]
  tx[,`:=`(ybar=ymax+.16*pad,ytext=ymax+.27*pad,
    label=paste0('Paired t: p=',fmt_p(p),'\nBH q=',fmt_p(q),' | n=',n_pairs))]
  p<-ggplot(boxes,aes(position,expression))+
   geom_boxplot(aes(group=position),width=STYLE$box_width,fill='#EFEFEF',color='#555555',linewidth=.4,outlier.shape=NA)+
   geom_segment(data=pairs,aes(x=1+offset,xend=2+offset,y=pre,yend=expression,color=response),inherit.aes=FALSE,linewidth=STYLE$line_width,alpha=STYLE$line_alpha)+
   geom_point(aes(x=position+offset,color=response,shape=response),size=STYLE$point_size,alpha=STYLE$point_alpha,show.legend=TRUE)+
   geom_segment(data=tx,aes(x=1,xend=2,y=ybar,yend=ybar),inherit.aes=FALSE,linewidth=.35)+
   geom_text(data=tx,aes(x=1.5,y=ytext,label=label),inherit.aes=FALSE,size=2.6,vjust=0,lineheight=1.15)+
   facet_wrap(~gene,nrow=1,scales='free_y')+
   scale_x_continuous(breaks=c(1,2),labels=c('Pre',tm),limits=c(.6,2.4))+
   scale_y_continuous(expand=expansion(mult=c(.06,.20)))+
   scale_color_manual(values=STYLE$response_colors,drop=FALSE,labels=c(R='Responder',NR='Non-responder'),limits=c('R','NR'))+
   scale_shape_manual(values=STYLE$point_shapes,drop=FALSE,labels=c(R='Responder',NR='Non-responder'),limits=c('R','NR'))+
   labs(x='Treatment node',y='Tumor pseudobulk log2(CPM + 1)',title=paste0(co,' | Pre vs ',tm),
    subtitle='Same patients on both sides; jittered points retain paired connections',
    caption=paste0('Box: median and IQR; whiskers: 1.5 IQR. p: paired t; q: cohort-wide BH correction.\n',if(nrow(pairs)/3==2)'Only 2 pairs: box summaries and t-test estimates are highly unstable.' else note))
  save_fig(p,paste0('Fig3_',co,'_Pre_vs_',tm,'_paired_boxplot'),STYLE$box_width_mm,STYLE$box_height_mm)
 }
}
write_json(STYLE,file.path(panel,'figure/style_used.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(panel,'figure/R_plot_sessionInfo.txt'))
