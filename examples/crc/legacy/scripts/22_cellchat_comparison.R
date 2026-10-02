# Compare exact CellChat pathway outputs using independent patients. Internal
# CellChat cell-label permutation P values gate pathway edges but are NOT used
# as high/low patient P values. Missing cell populations remain missing, not zero.
root<-normalizePath(file.path(dirname(sub('--file=','',grep('--file=',commandArgs(),value=TRUE)[1])),'..'),winslash='/')
.libPaths(c(file.path(root,'runtime/R-library'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(ggplot2)})
cfg<-fromJSON(file.path(root,'cellchat_config.json'));set.seed(cfg$seed)
jobs<-fread(file.path(root,'data/CellChat/jobs.csv'));hl<-fread(file.path(root,'tables/patient_high_low_groups.csv'))
db<-readRDS(file.path(root,'data/CellChat/CellChatDB_protein_locked.rds'));paths<-sort(c(unique(db$interaction$pathway_name),'MHC_I_classical_CD8'))
records<-list();lr<-list();audit<-list();networks<-list()
for(k in seq_len(nrow(jobs))){
 j<-jobs[k];f<-file.path(root,'derived/CellChat',paste0(j$job,'.rds'));stopifnot(file.exists(f));z<-readRDS(f)
 lev<-names(z$cell_counts)
 audit[[k]]<-data.table(job=j$job,study=j$study,patient=j$patient,n_cells=sum(z$cell_counts),n_groups=length(lev),n_LR_tested=dim(z$net$prob)[3],n_pathways_detected=length(z$netP$pathways))
 # Export all nonzero significant LR edges touching malignant cells for audit.
 ix<-which(z$net$prob>0&z$net$pval<=.05,arr.ind=TRUE)
 if(nrow(ix)){
  nm<-dimnames(z$net$prob);rr<-data.table(study=j$study,patient=j$patient,sample=j$sample,source=nm[[1]][ix[,1]],target=nm[[2]][ix[,2]],interaction=nm[[3]][ix[,3]],prob=z$net$prob[ix],cellchat_internal_p=z$net$pval[ix])
  lr[[k]]<-rr[source=='Malignant'|target=='Malignant']
 }
 for(a in lev)for(b in lev){
  if(a!='Malignant'&&b!='Malignant')next
  pp<-setNames(rep(0,length(paths)),paths)
  if(length(z$netP$pathways))pp[z$netP$pathways]<-z$netP$prob[a,b,]
  # Separate classical HLA-A/B/C--CD8A/B edges from the broad database MHC-I
  # pathway, which also contains nonclassical HLA and stress-ligand/NK edges.
  classical<-rownames(z$LR$LRsig)[z$LR$LRsig$ligand%in%c('HLA-A','HLA-B','HLA-C')&z$LR$LRsig$receptor%in%c('CD8A','CD8B')]
  classical<-intersect(classical,dimnames(z$net$prob)[[3]])
  if(length(classical))pp['MHC_I_classical_CD8']<-sum(z$net$prob[a,b,classical]*(z$net$pval[a,b,classical]<=.05))
  total<-sum(z$netP$prob)
  records[[length(records)+1]]<-data.table(study=j$study,patient=j$patient,sample=j$sample,source=a,target=b,pathway=paths,strength=as.numeric(pp),relative_strength=as.numeric(pp)/max(total,1e-15),n_source=as.integer(z$cell_counts[a]),n_target=as.integer(z$cell_counts[b]))
  networks[[length(networks)+1]]<-data.table(study=j$study,patient=j$patient,source=a,target=b,weight=z$net$weight[a,b])
 }
}
d<-rbindlist(records);fwrite(d,file.path(root,'tables/cellchat_sample_pathway_strength.csv'))
fwrite(rbindlist(lr),file.path(root,'tables/cellchat_significant_LR_edges.csv.gz'))
fwrite(rbindlist(audit),file.path(root,'audit/cellchat_completed_samples.csv'))
patient<-d[,.(strength=mean(strength),relative_strength=mean(relative_strength),n_samples=.N),by=.(study,patient,source,target,pathway)]
fwrite(patient,file.path(root,'tables/cellchat_patient_pathway_strength.csv'))
Ttypes<-c('CD8_T','CD4_T','Treg')
focus<-patient[(source=='Malignant'&target%in%Ttypes)|(target=='Malignant'&source%in%Ttypes)]
results<-list();elig<-list()
for(st in unique(focus$study))for(g in unique(hl$gene)){
 h<-hl[study==st&gene==g,.(patient,expression_group,target_mean)]
 u<-merge(focus[study==st],h,by='patient')
 for(edge in unique(paste(u$source,u$target,sep='|'))){
  ed<-strsplit(edge,'|',fixed=TRUE)[[1]];v<-u[source==ed[1]&target==ed[2]]
  pts<-sort(unique(v$patient));hh<-h[match(pts,patient)];isH<-hh$expression_group=='High';nh<-sum(isH);nl<-sum(!isH)
  elig[[length(elig)+1]]<-data.table(study=st,gene=g,source=ed[1],target=ed[2],n_high=nh,n_low=nl,eligible=min(nh,nl)>=cfg$min_patients_each_high_low)
  if(min(nh,nl)<cfg$min_patients_each_high_low)next
  # Reuse patient-label permutations across pathways on the identical patient set.
  exact_permutation<-choose(length(pts),nh)<=cfg$permutations
  if(exact_permutation){
   combos<-combn(seq_along(pts),nh);W<-matrix(-1/nl,nrow=length(pts),ncol=ncol(combos))
   for(ii in seq_len(ncol(combos)))W[combos[,ii],ii]<-1/nh
  }else W<-replicate(cfg$permutations,{b<-sample(isH);as.numeric(b)/nh-as.numeric(!b)/nl})
  w<-as.numeric(isH)/nh-as.numeric(!isH)/nl
  for(metric in c('strength','relative_strength')){
   wide<-dcast(v,pathway~patient,value.var=metric)
   mat<-as.matrix(wide[,..pts]);rownames(mat)<-wide$pathway
   ranks<-t(apply(mat,1,rank,ties.method='average'));obs<-as.numeric(ranks%*%w);null<-ranks%*%W
   permP<-if(exact_permutation)rowMeans(abs(null)>=abs(obs)-1e-10) else (1+rowSums(abs(null)>=abs(obs)-1e-10))/(cfg$permutations+1)
   ct<-lapply(seq_len(nrow(mat)),function(i){x<-mat[i,];if(sd(x)==0)c(rho=NA_real_,p=1) else{exact_rho<-length(x)<=9&&!anyDuplicated(x)&&!anyDuplicated(hh$target_mean);a<-cor.test(x,hh$target_mean,method='spearman',exact=exact_rho);c(rho=unname(a$estimate),p=a$p.value)}})
   corrs<-do.call(rbind,ct)
   results[[length(results)+1]]<-data.table(study=st,gene=g,source=ed[1],target=ed[2],pathway=rownames(mat),metric=metric,n_high=nh,n_low=nl,permutation_mode=if(exact_permutation)'exact' else 'Monte_Carlo',permutations_used=ncol(W),
       mean_high=rowMeans(mat[,isH,drop=FALSE]),mean_low=rowMeans(mat[,!isH,drop=FALSE]),
       median_high=apply(mat[,isH,drop=FALSE],1,median),median_low=apply(mat[,!isH,drop=FALSE],1,median),
       rank_mean_difference=obs,p_patient_permutation=permP,rho=corrs[,'rho'],p_spearman=corrs[,'p'])
  }
 }
}
r<-rbindlist(results);r[,q_patient_permutation:=p.adjust(p_patient_permutation,'BH'),by=.(study,metric)]
r[,q_spearman:=p.adjust(p_spearman,'BH'),by=.(study,metric)]
r[,`:=`(q_focused_MHC_permutation=NA_real_,q_focused_MHC_spearman=NA_real_)]
r[pathway=='MHC_I_classical_CD8'&source=='Malignant'&target=='CD8_T',q_focused_MHC_permutation:=p.adjust(p_patient_permutation,'BH'),by=.(study,metric)]
r[pathway=='MHC_I_classical_CD8'&source=='Malignant'&target=='CD8_T',q_focused_MHC_spearman:=p.adjust(p_spearman,'BH'),by=.(study,metric)]
fwrite(r,file.path(root,'tables/cellchat_patient_high_low_tests.csv'));fwrite(rbindlist(elig),file.path(root,'tables/cellchat_comparison_eligibility.csv'))
# Plot a fixed immunology panel, with no outcome-dependent pathway selection.
panel<-intersect(c('MHC_I_classical_CD8','MHC-I','MHC-II','CD80','CD86','CD28','CD45','PD-L1','PDL2','GALECTIN','NECTIN','PVR','CXCL','CCL','IFN-II','TNF','TGFb','IL2','IL10'),paths)
theme_set(theme_classic(base_size=9)+theme(plot.title=element_text(size=11),legend.position='bottom'))
saveplot<-function(p,name,w=183,h=170){ggsave(file.path(root,'figures',paste0(name,'.pdf')),p,width=w,height=h,units='mm',device=cairo_pdf);ggsave(file.path(root,'figures',paste0(name,'.png')),p,width=w,height=h,units='mm',dpi=300,bg='white')}
for(st in unique(r$study)){
 v<-r[study==st&metric=='strength'&pathway%in%panel];v[,edge:=paste(source,'to',target)]
 v[,label:=ifelse(is.na(rho),'-',sprintf('%.2f%s',rho,ifelse(q_spearman<.05,'*','')))]
 p<-ggplot(v,aes(gene,pathway,fill=rho))+geom_tile(color='white')+geom_text(aes(label=label),size=2.1)+facet_wrap(~edge,ncol=3)+scale_fill_gradient2(low='#4C78A8',mid='white',high='#C76D5D',limits=c(-1,1),na.value='grey95')+labs(title=paste('CellChat: tumor expression and T-cell communication,',st),subtitle='Patient-level rho; * exploratory BH q < 0.05; dash = constant score, rho not estimable',x=NULL,y=NULL)
 saveplot(p,paste0('09_CellChat_pathway_',st),h=205)
}
# Raw patient points for predefined tumor-to-CD8 MHC-I signaling, if in DB.
# Each patient intentionally has two pathway rows and three target-gene labels.
v<-merge(patient[source=='Malignant'&target=='CD8_T'&pathway%in%c('MHC-I','MHC_I_classical_CD8')],hl,by=c('study','patient'),allow.cartesian=TRUE)
if(nrow(v)){
 p<-ggplot(v,aes(expression_group,strength,color=expression_group))+geom_boxplot(outlier.shape=NA,width=.5)+geom_point(position=position_jitter(width=.10,height=0,seed=cfg$seed),size=1.5)+facet_grid(study+pathway~gene,scales='free_y')+scale_color_manual(values=c(High='#CC7956',Low='#287F90'))+labs(title='CellChat MHC-I: malignant cells to CD8 T cells',subtitle='One point per patient; computational communication score, not antigen-specific presentation',x='Tumor target expression group',y='CellChat pathway strength')
 saveplot(p,'10_CellChat_MHC_I_patient_points',h=200)
}
nw<-merge(rbindlist(networks),hl,by=c('study','patient'),allow.cartesian=TRUE);nw<-nw[,.(mean_weight=mean(weight)),by=.(study,gene,expression_group,source,target)]
fwrite(nw,file.path(root,'tables/cellchat_high_low_network_means.csv'))
for(st in unique(nw$study)){
 v<-nw[study==st];p<-ggplot(v,aes(source,target,fill=log1p(mean_weight)))+geom_tile(color='white')+facet_grid(expression_group~gene)+scale_fill_gradient(low='white',high='#287F90')+theme(axis.text.x=element_text(angle=55,hjust=1,size=6),axis.text.y=element_text(size=6))+labs(title=paste('CellChat tumor-centered network:',st),subtitle='Group means across patients; only malignant-cell incoming/outgoing edges shown',x='Sender',y='Receiver',fill='log1p(mean weight)')
 saveplot(p,paste0('11_CellChat_network_',st),h=160)
}
cat('COMPARISON COMPLETE',nrow(jobs),'samples;',nrow(r),'tests (including relative-score sensitivity)\n')
print(head(r[metric=='strength'&q_patient_permutation<.05][order(q_patient_permutation)],15))
