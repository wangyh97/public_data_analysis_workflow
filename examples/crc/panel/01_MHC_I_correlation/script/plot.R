# Plotting only: reads source/*.csv; never reads raw matrices or redefines tests.
# Run: Rscript plot.R. All graphics are made with R/ggplot2 and editable PDF text.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite)})
script<-sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])
panel_dir<-normalizePath(file.path(dirname(script),'..'),winslash='/')

# ================= USER STYLE INTERFACE / 自定义样式 =================
STYLE<-list(
  font='Arial',font_size=9,title_size=11,strip_size=9,
  negative='#326A9F',zero='#FAFAFA',positive='#B64B45',missing='#D9D9D9',
  cohort_colors=c(GSE236581='#326A9F',GSE205506='#B86B38'),
  heat_width_mm=183,heat_height_mm=108,heat_number_size=3.1,
  show_rho=TRUE,show_q_stars=TRUE,star_cutoffs=c(.05,.01,.001),
  # Cell text can instead be 'q' or 'n'; stars always refer to BH q, never raw p.
  heat_text='rho',rho_limits=c(-1,1),
  scatter_width_mm=183,scatter_height_mm=88,point_size=2.2,point_alpha=.85,
  point_shape=21,point_stroke=.45,label_patients=FALSE,
  trend_line=TRUE,trend_method='lm',trend_se=FALSE,trend_width=.5,
  # lm is a descriptive guide on expression axes; statistical test is Spearman.
  dpi=600,formats=c('pdf','png'),
  target_order=c('METTL3','STRAP','PTBP1'),
  cohort_order=c('GSE236581','GSE205506'),
  module_order=c('MHC_I_structure','Immunoproteasome','Peptide_transport',
    'ER_peptide_trimming','Folding_and_loading','Transcriptional_regulators'),
  module_labels=c(MHC_I_structure='MHC-I structure',Immunoproteasome='Immunoproteasome',
    Peptide_transport='Peptide transport into ER',ER_peptide_trimming='ER peptide trimming',
    Folding_and_loading='Folding and peptide loading',Transcriptional_regulators='Transcriptional regulators'),
  structural_order=c('HLA-A','HLA-B','HLA-C','B2M')
)
# Change values above, rerun this script; no processing rerun is needed for style.
# ================= END STYLE INTERFACE ==============================
dir.create(file.path(panel_dir,'figure'),showWarnings=FALSE)
theme_set(theme_classic(base_size=STYLE$font_size,base_family=STYLE$font)+
  theme(plot.title=element_text(size=STYLE$title_size,face='bold',margin=margin(b=6)),
    plot.subtitle=element_text(size=8,margin=margin(b=9)),
    plot.caption=element_text(size=7,hjust=0,margin=margin(t=9)),
    strip.background=element_blank(),strip.text=element_text(size=STYLE$strip_size,face='bold'),
    axis.text=element_text(color='#303030'),axis.line=element_line(linewidth=.35),
    axis.ticks=element_line(linewidth=.35),panel.spacing=grid::unit(8,'mm'),
    plot.margin=margin(9,11,9,9),legend.position='bottom'))
save_figure<-function(p,name,w,h){
 for(ext in STYLE$formats){
   dest<-file.path(panel_dir,'figure',paste0(name,'.',ext))
   if(ext=='pdf'){
     grDevices::cairo_pdf(dest,width=w/25.4,height=h/25.4,family=STYLE$font);print(p);dev.off()
   }
   # Optional .svg export needs svglite; default .pdf / .png need base R only.
   else if(ext=='svg'){
     if(!requireNamespace('svglite',quietly=TRUE))stop('Optional SVG export requires svglite')
     svglite::svglite(dest,width=w/25.4,height=h/25.4);print(p);dev.off()
   }
   else if(ext=='png')ggsave(dest,p,width=w,height=h,units='mm',dpi=STYLE$dpi,device=grDevices::png,type='cairo')
   else if(ext=='tiff')ggsave(dest,p,width=w,height=h,units='mm',dpi=STYLE$dpi,compression='lzw',type='cairo')
   else stop('Supported formats: pdf, png, tiff, svg')
 }
}
qstars<-function(q)ifelse(is.na(q),'',ifelse(q<STYLE$star_cutoffs[3],'***',ifelse(q<STYLE$star_cutoffs[2],'**',ifelse(q<STYLE$star_cutoffs[1],'*',''))))
fmt<-function(x)ifelse(is.na(x),'NA',formatC(x,format='g',digits=2))
stats<-fread(file.path(panel_dir,'source/correlation_statistics.csv'))
pts<-fread(file.path(panel_dir,'source/patient_values.csv'))
stats[,cohort:=factor(cohort,levels=STYLE$cohort_order)]
stats[,target:=factor(target,levels=STYLE$target_order)]
h<-stats[family=='MHC_modules']
h[,outcome:=factor(outcome,levels=rev(STYLE$module_order))]
h[,label:=if(STYLE$heat_text=='q')paste0('q=',fmt(q)) else if(STYLE$heat_text=='n')paste0('n=',n) else if(STYLE$show_rho)sprintf('%.2f',rho) else '']
if(STYLE$show_q_stars)h[,label:=paste0(label,qstars(q))]
h[is.na(rho),label:='NA']
ns<-h[,.(n_label=paste(sort(unique(n)),collapse='-')),by=cohort]
cohort_labels<-setNames(paste0(ns$cohort,'  |  n = ',ns$n_label),ns$cohort)
p<-ggplot(h,aes(target,outcome,fill=rho))+geom_tile(color='white',linewidth=1.1)+
 geom_text(aes(label=label),size=STYLE$heat_number_size,color='#202020')+
 facet_grid(~cohort,labeller=labeller(cohort=cohort_labels))+
 scale_fill_gradient2(low=STYLE$negative,mid=STYLE$zero,high=STYLE$positive,midpoint=0,limits=STYLE$rho_limits,na.value=STYLE$missing,breaks=c(-1,-.5,0,.5,1))+
 scale_y_discrete(labels=STYLE$module_labels,expand=c(0,0))+scale_x_discrete(expand=c(0,0),position='top')+
 labs(x=NULL,y=NULL,fill='Spearman rho',title='Tumor gene expression and MHC-I programs',
 subtitle='Pretreatment CRC tumors  |  Patient-level correlations',
 caption='Scores and target expression: malignant-labelled cells. * q < 0.05; ** q < 0.01; *** q < 0.001.\nBH correction: 18 module tests per cohort. GSE205506 uses scCT-DB malignant annotations.')+
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(face='italic'),
   legend.key.width=grid::unit(13,'mm'),legend.key.height=grid::unit(3,'mm'))
save_figure(p,'Fig1_MHC_I_module_heatmap',STYLE$heat_width_mm,STYLE$heat_height_mm)

# One independent file per target/structural-gene pair, with both cohorts.
for(g in STYLE$target_order)for(y_gene in STYLE$structural_order){
 d<-pts[target==g & outcome==y_gene]
 st<-stats[target==g & outcome==y_gene]
 st[,facet_label:=paste0(cohort,'  |  n = ',n,'\nrho = ',sprintf('%.2f',rho),'   p = ',fmt(p),'   q = ',fmt(q))]
 d[,facet_label:=factor(st$facet_label[match(cohort,st$cohort)],levels=st$facet_label[order(st$cohort)])]
 p<-ggplot(d,aes(x,y))+geom_point(aes(fill=cohort),shape=STYLE$point_shape,size=STYLE$point_size,
   alpha=STYLE$point_alpha,stroke=STYLE$point_stroke,color='white')+
   facet_wrap(~facet_label,nrow=1,scales='free')+
   scale_fill_manual(values=STYLE$cohort_colors,guide='none')+
   scale_x_continuous(expand=expansion(mult=c(.08,.08)))+scale_y_continuous(expand=expansion(mult=c(.08,.1)))+
   labs(x=paste0(g,' expression'),y=paste0(y_gene,' expression'),
     title=paste0(g,' and ',y_gene),
     caption='Each point is one patient; axes show malignant-cell mean log1p(CP10K).\nLine: descriptive linear fit. Spearman test; BH q across 12 structural-gene tests per cohort.')
 if(STYLE$trend_line)p<-p+geom_smooth(method=STYLE$trend_method,formula=y~x,se=STYLE$trend_se,color='#666666',linewidth=STYLE$trend_width)
 if(STYLE$label_patients)p<-p+geom_text(aes(label=patient),size=2.4,vjust=-.7,check_overlap=FALSE)
 save_figure(p,paste0('Fig1_scatter_',g,'_',y_gene),STYLE$scatter_width_mm,STYLE$scatter_height_mm)
}
writeLines(capture.output(sessionInfo()),file.path(panel_dir,'figure/R_plot_sessionInfo.txt'))
write_json(STYLE,file.path(panel_dir,'figure/style_used.json'),pretty=TRUE,auto_unbox=TRUE)
