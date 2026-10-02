# R-only plotting from portable patient-level source CSV. No raw data required.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite)})
script<-sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])
panel_dir<-normalizePath(file.path(dirname(script),'..'),winslash='/')
# ================= USER STYLE INTERFACE / 自定义样式 =================
STYLE<-list(font='Arial',font_size=9,title_size=11,strip_size=9,
 negative='#326A9F',zero='#FAFAFA',positive='#B64B45',missing='#D9D9D9',
 width_mm=183,height_mm=122,dpi=600,formats=c('pdf','png'),
 number_size=3.1,show_q_stars=TRUE,heat_text='rho',rho_limits=c(-1,1),
 target_order=c('METTL3','STRAP','PTBP1'),cohort_order=c('GSE236581','GSE205506'),
 outcome_order=c('CD8_T_Cytotoxic','CD8_T_Dysfunction','CD8_T_Stem_memory',
   'All_T_Cytotoxic','All_T_Dysfunction','All_T_Stem_memory','All_T_Treg'),
 outcome_labels=c(CD8_T_Cytotoxic='CD8 T: cytotoxicity',CD8_T_Dysfunction='CD8 T: dysfunction',
   CD8_T_Stem_memory='CD8 T: stem / memory',All_T_Cytotoxic='All T: cytotoxicity',
   All_T_Dysfunction='All T: dysfunction',All_T_Stem_memory='All T: stem / memory',
   All_T_Treg='All T: Treg program'))
# Edit the values above to customize font, size, palette, order and text.
# heat_text supports 'rho', 'q', or 'n'; stars always indicate BH q.
# ================= END STYLE INTERFACE ==============================
theme_set(theme_classic(base_size=STYLE$font_size,base_family=STYLE$font)+theme(
 plot.title=element_text(size=STYLE$title_size,face='bold',margin=margin(b=6)),
 plot.subtitle=element_text(size=8,margin=margin(b=9)),
 plot.caption=element_text(size=7,hjust=0,margin=margin(t=9)),
 strip.background=element_blank(),strip.text=element_text(size=STYLE$strip_size,face='bold'),
 axis.text=element_text(color='#303030'),axis.line=element_blank(),axis.ticks=element_blank(),
 panel.spacing=grid::unit(8,'mm'),plot.margin=margin(9,11,9,9),legend.position='bottom',
 legend.key.width=grid::unit(13,'mm'),legend.key.height=grid::unit(3,'mm')))
d<-fread(file.path(panel_dir,'source/correlation_statistics.csv'))
d[,cohort:=factor(cohort,levels=STYLE$cohort_order)]
d[,target:=factor(target,levels=STYLE$target_order)]
d[,outcome:=factor(outcome,levels=rev(STYLE$outcome_order))]
d[,label:=if(STYLE$heat_text=='q')paste0('q=',formatC(q,digits=2,format='g')) else if(STYLE$heat_text=='n')paste0('n=',n) else sprintf('%.2f',rho)]
if(STYLE$show_q_stars)d[,label:=paste0(label,ifelse(q<.001,'***',ifelse(q<.01,'**',ifelse(q<.05,'*',''))))]
d[is.na(rho),label:='NA']
ns<-d[,.(n_label=paste(sort(unique(n)),collapse='-')),by=cohort]
cohort_labels<-setNames(paste0(ns$cohort,'  |  n = ',ns$n_label),ns$cohort)
p<-ggplot(d,aes(target,outcome,fill=rho))+geom_tile(color='white',linewidth=1.1)+
 geom_text(aes(label=label),size=STYLE$number_size,color='#202020')+
 # Separator between CD8-only and all-T programs; cells are overlapping sets.
 geom_hline(yintercept=4.5,color='white',linewidth=2)+
 facet_grid(~cohort,labeller=labeller(cohort=cohort_labels))+
 scale_fill_gradient2(low=STYLE$negative,mid=STYLE$zero,high=STYLE$positive,midpoint=0,limits=STYLE$rho_limits,na.value=STYLE$missing,breaks=c(-1,-.5,0,.5,1))+
 scale_x_discrete(position='top',expand=c(0,0))+scale_y_discrete(labels=STYLE$outcome_labels,expand=c(0,0))+
 labs(x=NULL,y=NULL,fill='Spearman rho',title='Tumor gene expression and T-cell states',
 subtitle='Pretreatment CRC tumors  |  Matched within the same tumor specimen',
 caption='Target expression: malignant-labelled cells; state scores: indicated T-cell compartment.\n* q < 0.05; ** q < 0.01; *** q < 0.001. BH correction: 21 tests per cohort.\nAll-T and CD8-T scores overlap. GSE205506 uses scCT-DB malignant annotations.')+
 theme(axis.text.x=element_text(face='italic'))
for(ext in STYLE$formats){
 dest<-file.path(panel_dir,'figure',paste0('Fig2_T_cell_state_heatmap.',ext))
 if(ext=='pdf'){
   grDevices::cairo_pdf(dest,width=STYLE$width_mm/25.4,height=STYLE$height_mm/25.4,family=STYLE$font);print(p);dev.off()
 }
 # Optional .svg export; default .pdf and .png use base R graphics devices.
 else if(ext=='svg'){
   if(!requireNamespace('svglite',quietly=TRUE))stop('Optional SVG export requires svglite')
   svglite::svglite(dest,width=STYLE$width_mm/25.4,height=STYLE$height_mm/25.4);print(p);dev.off()
 }
 else if(ext=='png')ggsave(dest,p,width=STYLE$width_mm,height=STYLE$height_mm,units='mm',dpi=STYLE$dpi,device=grDevices::png,type='cairo')
 else if(ext=='tiff')ggsave(dest,p,width=STYLE$width_mm,height=STYLE$height_mm,units='mm',dpi=STYLE$dpi,compression='lzw',type='cairo')
 else stop('Supported formats: pdf, png, tiff, svg')
}
write_json(STYLE,file.path(panel_dir,'figure/style_used.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(panel_dir,'figure/R_plot_sessionInfo.txt'))
