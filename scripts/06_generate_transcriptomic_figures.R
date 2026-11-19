required_packages <- c('data.table','ggplot2','patchwork','ggrepel','ggnewscale','ragg','svglite','scales')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_packages)) stop('Install required R packages: ', paste(missing_packages, collapse=', '))

suppressPackageStartupMessages({
 library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
 library(ggnewscale); library(ragg); library(svglite); library(scales)
})

data_dir<-file.path(getwd(),'results/plot_data')
out<-file.path(getwd(),'figures')
for(d in c('transcriptomic_PCA','pathway_enrichment','SREBP1_lipogenic_signaling','source_data'))dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
models <- c('OVCAR5','OVCAR8','ID8')
cellcols <- c(OVCAR5='#E9B4B4',OVCAR8='#B8D7F0',ID8='#B9DDBE')
treatcols <- c(Vehicle='#8F969C',NB25='#E76F51')
linecols <- c(OVCAR5='#AF6672',OVCAR8='#477EAC',ID8='#56936E')
blue <- '#3D70A2'; orange <- '#D97158'
savefig <- function(p,folder,stem,w,h) {
 stem<-sub('^Figure_5A_design$','experimental_design',stem)
 stem<-sub('^Figure_5[B-D]_PCA_','PCA_',stem)
 stem<-sub('^Figure_5E_selected_Hallmark$','selected_Hallmark_enrichment',stem)
 stem<-sub('^Figure_5_final$','transcriptomic_profiling',stem)
 stem<-sub('^Figure_6[A-C]_module_GSEA_','SREBP1_module_GSEA_',stem)
 stem<-sub('^Figure_6D_lipogenic_replicates$','lipogenic_gene_replicates',stem)
 stem<-sub('^Figure_6E_five_program_comparison$','lipid_adipogenesis_programs',stem)
 stem<-sub('^Figure_6F_sterol_gene_effects$','sterol_gene_effects',stem)
 stem<-sub('^Figure_6_final$','SREBP1_lipogenic_signaling',stem)
 folder<-if(folder=='Figure_5'){if(grepl('Hallmark',stem))'pathway_enrichment' else 'transcriptomic_PCA'}else 'SREBP1_lipogenic_signaling'
 b<-file.path(out,folder,stem)
 ggsave(paste0(b,'.pdf'),p,width=w,height=h,device=grDevices::cairo_pdf,bg='white',limitsize=FALSE)
 ggsave(paste0(b,'.svg'),p,width=w,height=h,device=svglite::svglite,bg='white',limitsize=FALSE)
 ggsave(paste0(b,'.png'),p,width=w,height=h,dpi=400,device=ragg::agg_png,bg='white',limitsize=FALSE)
}
stars <- function(q) fifelse(is.na(q),'',fifelse(q<.001,'***',fifelse(q<.01,'**',fifelse(q<.05,'*',''))))
base_theme <- theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(size=16,face='bold',hjust=0),plot.subtitle=element_text(size=11),plot.title.position='plot',plot.margin=margin(5,5,5,5))
meta <- fread(file.path(data_dir,'sample_metadata_confirmed_pairs.csv'))
meta[,cell_line:=factor(cell_line,levels=models)]
meta[,treatment:=factor(treatment,levels=c('CON','NB25'))]
stopifnot(nrow(meta)==18, all(meta[,table(cell_line,treatment)]==3), uniqueN(meta$sample_id)==18)
fwrite(meta[,.(sample_id,cell_line,treatment,replicate,batch)],file.path(out,'source_data','all_18_sample_metadata.csv'))
# Figure 5A: sample design from confirmed metadata.
design <- meta[,.(n=.N),by=.(cell_line,treatment)]
design[,treatment_label:=factor(ifelse(treatment=='CON','Vehicle','NB25'),levels=c('Vehicle','NB25'))]
design[,cell_line:=factor(cell_line,levels=rev(models))]
p5A <- ggplot(design,aes(treatment_label,cell_line,fill=treatment_label))+
 geom_tile(width=.84,height=.72,color='white',linewidth=1)+
 geom_text(aes(label=paste0('n = ',n)),size=12/.pt,fontface='bold')+
 scale_fill_manual(values=treatcols)+scale_x_discrete(labels=c(Vehicle='Vehicle',NB25='NDNB25'))+
 labs(title='A. RNA-seq experimental design',subtitle='OVCAR5 and OVCAR8 (human), ID8 (mouse); 24 h NDNB25 or vehicle; biological replicates',x=NULL,y=NULL)+
 base_theme+theme(panel.grid=element_blank(),legend.position='none',axis.text=element_text(size=11,face='bold'))
savefig(p5A,'Figure_5','Figure_5A_design',7.5,1.8)
fwrite(design,file.path(out,'source_data','experimental_design.csv'))
# Figure 5B-D: exact archived within-cell-line VST PCA coordinates.
vv <- fread(file.path(data_dir,'VST_PCA_variance.csv'))
pca <- lapply(seq_along(models),function(i){
 z <- models[i]
 d <- fread(file.path(data_dir,paste0(z,'_paired_VST_PCA_scores.csv')))
 d[,cell_line:=z]; d[,replicate:=as.integer(sub('.*_','',label))]
 d <- merge(d,vv,by='cell_line')
 fwrite(d,file.path(out,'source_data',paste0(z,'_PCA.csv')))
 d[,treatment:=factor(treatment,levels=c('CON','NB25'))]
 p <- ggplot(d,aes(PC1,PC2,color=treatment,shape=treatment))+
  geom_point(size=3.5,stroke=.55)+
  ggrepel::geom_text_repel(aes(label=replicate),size=10/.pt,color='#24282A',seed=20260924,min.segment.length=0,point.padding=.12,box.padding=.1,segment.color='#9DA6AD',segment.size=.2,max.overlaps=Inf,show.legend=FALSE)+
  scale_color_manual(values=c(CON=treatcols[['Vehicle']],NB25=treatcols[['NB25']]),labels=c(CON='Vehicle',NB25='NDNB25'))+
  scale_shape_manual(values=c(CON=16,NB25=17),labels=c(CON='Vehicle',NB25='NDNB25'))+
  scale_x_continuous(expand=expansion(mult=c(.12,.14)))+scale_y_continuous(expand=expansion(mult=c(.12,.14)))+
  labs(title=paste0(LETTERS[i+1],'. ',z),x=sprintf('PC1 (%.1f%%)',d$PC1_percent[1]),y=sprintf('PC2 (%.1f%%)',d$PC2_percent[1]))+
  base_theme+theme(plot.title=element_text(size=14,face='bold'),axis.title=element_text(size=11,margin=margin(2,2,2,2)),axis.text=element_text(size=10),
   panel.grid.major=element_line(color='#E8ECF0',linewidth=.25),aspect.ratio=1,legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=10),plot.margin=margin(3,3,3,3))
 savefig(p,'Figure_5',paste0('Figure_5',LETTERS[i+1],'_PCA_',z),3.3,3.45)
 p
})
pca_row <- wrap_plots(pca,nrow=1,guides='collect') & theme(legend.position='bottom')
# Figure 5E: explicitly selected after exploration; not called prespecified.
hall <- rbindlist(lapply(models,function(z){
 d<-fread(file.path(data_dir,paste0(z,'_Hallmark_only_seeded_multilevel_GSEA.csv')))
 d[,cell_line:=z]; d
}))
sel <- c('HALLMARK_CHOLESTEROL_HOMEOSTASIS','HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_DNA_REPAIR','HALLMARK_APOPTOSIS','HALLMARK_E2F_TARGETS',
 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_INFLAMMATORY_RESPONSE',
 'HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_INTERFERON_ALPHA_RESPONSE')
labels <- c('Cholesterol homeostasis','Fatty-acid metabolism','Oxidative phosphorylation',
 'Unfolded-protein response','DNA repair','Apoptosis','E2F targets','EMT',
 'TNF-alpha/NF-kappaB signaling','Inflammatory response','IL-6/JAK/STAT3 signaling','Interferon-alpha response')
yy <- c(14,13,12,10.5,9.5,8.5,7.5,6.5,5,4,3,2)
g5 <- hall[pathway %in% sel]
stopifnot(nrow(g5)==36)
g5[,display:=labels[match(pathway,sel)]]
g5[,group:=c(rep('Metabolism',3),rep('Stress / cell state',5),rep('Inflammatory signaling',4))[match(pathway,sel)]]
g5[,y:=yy[match(pathway,sel)]]
g5[,cell_line:=factor(cell_line,levels=models)]
g5[,txt:=paste0(sprintf('%.2f',NES),stars(padj))]
fwrite(g5[,.(pathway,display,group,cell_line,NES,nominal_p=pval,FDR=padj,txt,y,size)],file.path(out,'source_data','selected_Hallmark.csv'))
p5E <- ggplot(g5,aes(cell_line,y,fill=NES))+geom_tile(height=.92,width=.98,color='white',linewidth=.5)+
 geom_text(aes(label=txt),size=11/.pt,fontface='bold')+
 scale_fill_gradient2(low=blue,mid='white',high=orange,midpoint=0,limits=c(-3,3),name='NES')+
 scale_x_discrete(position='top')+scale_y_continuous(breaks=yy,labels=labels,limits=c(1.5,14.5),expand=c(0,0))+
 labs(title='E. Selected Hallmark pathway enrichment',subtitle='NES shown; * FDR < 0.05, ** < 0.01, *** < 0.001. Pathways selected after exploratory review.',x=NULL,y=NULL)+
 base_theme+theme(panel.grid=element_blank(),axis.text.x=element_text(size=14,face='bold'),axis.text.y=element_text(size=10.5,face='bold'),
 legend.position='right',legend.title=element_text(size=10),legend.text=element_text(size=10),plot.margin=margin(4,4,4,3))
savefig(p5E,'Figure_5','Figure_5E_selected_Hallmark',7.5,5.5)
fig5 <- wrap_elements(full=p5A)/wrap_elements(full=pca_row)/wrap_elements(full=p5E)+plot_layout(heights=c(1.8,3.8,5.5))
savefig(fig5,'Figure_5','Figure_5_final',7.5,11.1)
# Figure 6A-C: running-score curves from archived ranks and audited fixed module, without rerunning GSEA.
mod <- fread(file.path(data_dir,'SREBP1_lipogenic_module_statistics.csv'))
mod[,model_order:=match(cell_line,models)]; setorder(mod,model_order); mod[,model_order:=NULL]
stopifnot(nrow(mod)==3,all(mod$NES<0),all(mod$BH_q_across_three_prespecified_tests<.05))
fwrite(mod,file.path(out,'source_data','SREBP1_module_statistics.csv'))
curves <- lapply(seq_along(models),function(i){
 z <- models[i]; q <- mod[cell_line==z]
 r <- fread(file.path(data_dir,paste0(z,'_ranked_GSEA.csv')))
 r <- r[!is.na(gene_label)&gene_label!=''&!is.na(stat)]
 r[,key:=toupper(gene_label)]
 r <- r[!duplicated(key)]
 setorder(r,-stat)
 genes <- toupper(strsplit(q$present_genes,';',fixed=TRUE)[[1]])
 hit <- r$key %in% genes
 stopifnot(sum(hit)==q$gene_set_size)
 wh <- abs(r$stat[hit]); weights <- numeric(nrow(r)); weights[hit] <- wh/sum(wh); weights[!hit] <- -1/(nrow(r)-sum(hit))
 running <- cumsum(weights)
 stopifnot(abs(min(running)-q$ES)<.015)
 curve <- data.table(cell_line=z,rank=seq_len(nrow(r)),running_ES=running,hit=hit,gene=r$key)
 fwrite(curve,file.path(out,'source_data',paste0(z,'_running_enrichment.csv')))
 p_curve <- ggplot(curve,aes(rank,running_ES))+geom_hline(yintercept=0,color='#BEC3C6',linewidth=.3)+
  geom_line(color='#318A68',linewidth=.85)+
  coord_cartesian(ylim=c(-.86,.24),xlim=c(1,nrow(r)),expand=FALSE)+
  labs(title=paste0(LETTERS[i],'. ',z),subtitle=sprintf('NES = %.2f\nBH q = %.4f',q$NES,q$BH_q_across_three_prespecified_tests),
   x=NULL,y='Running enrichment score')+
  theme_classic(base_size=11)+theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=11),
   axis.title.y=element_text(size=11),axis.text.y=element_text(size=10),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
   plot.margin=margin(4,5,0,4))
 rug <- curve[hit==TRUE,.(rank,gene)]
 stopifnot(nrow(rug)==q$gene_set_size,all(rug$gene %in% genes))
 p_rug <- ggplot(rug)+geom_segment(aes(x=rank,xend=rank,y=0,yend=1),color='#60676B',linewidth=.3)+
  coord_cartesian(xlim=c(1,nrow(r)),ylim=c(-.02,1.02),expand=FALSE)+
  scale_x_continuous(limits=c(1,nrow(r)),expand=c(0,0))+
  labs(x='Gene rank',y=NULL)+
  theme_classic(base_size=11)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),
   axis.title.x=element_text(size=11),axis.text.x=element_text(size=10),
   plot.margin=margin(0,5,4,4))
 p <- p_curve/p_rug+plot_layout(heights=c(3,.43))
 savefig(p,'Figure_6',paste0('Figure_6',LETTERS[i],'_module_GSEA_',z),3.3,3.3)
 p
})
p6AC <- wrap_plots(curves,nrow=1)+plot_annotation(title='A-C. SREBF1-associated lipogenic program enrichment',
 subtitle='Negative NES: module genes enrich toward NDNB25-decreased ranks; rug ticks mark gene hits.',
 theme=theme(plot.title=element_text(size=16,face='bold'),plot.subtitle=element_text(size=11),plot.title.position='plot',plot.margin=margin(3,3,3,3)))
# Figure 6D: eight lipogenic rows, separate Scd1 and Scd2; all 18 individual samples.
hm <- fread(file.path(data_dir,'selected_gene_VST_by_sample.csv'))
hm[,display_gene:=fifelse(display_gene=='SCD','SCD / Scd1',
 fifelse(display_gene=='SCD (Scd2)','Scd2 (mouse only)',display_gene))]
lipgenes <- c('ACLY','ACACA','FASN','SCD / Scd1','Scd2 (mouse only)','ELOVL6','FADS1','FADS2')
hm <- hm[display_gene %in% lipgenes,.(sample_id,display_gene,cell_line,z,vst,gene_id)]
m0 <- meta[,.(sample_id,cell_line,treatment,replicate)]
m0[,x:=match(as.character(cell_line),models)*6-6+replicate+ifelse(treatment=='NB25',3,0)]
grid <- rbindlist(lapply(seq_len(nrow(m0)),function(i)data.table(sample_id=m0$sample_id[i],display_gene=lipgenes)))
hm <- merge(grid,hm,by=c('sample_id','display_gene'),all.x=TRUE,sort=FALSE)
hm <- merge(hm,m0,by='sample_id',all.x=TRUE,suffixes=c('','_meta'))
stopifnot(nrow(hm)==18*length(lipgenes),uniqueN(hm$sample_id)==18)
hm[,y:=length(lipgenes)-match(display_gene,lipgenes)+1]
hm[,structural_NA:=display_gene=='Scd2 (mouse only)' & as.character(cell_line_meta) %in% c('OVCAR5','OVCAR8')]
hm[,value_status:=fifelse(structural_NA,'species_not_applicable',fifelse(is.na(z),'missing_processed_value','observed'))]
stopifnot(sum(hm$structural_NA)==12,all(!is.na(hm$z[!hm$structural_NA])))
fwrite(hm,file.path(out,'source_data','lipogenic_replicate_VST_z.csv'))
bars <- rbindlist(list(m0[,.(x,y=length(lipgenes)+2.05,category=as.character(cell_line))],
 m0[,.(x,y=length(lipgenes)+1.2,category=ifelse(treatment=='CON','Vehicle','NB25'))]))
bars[,category:=factor(category,levels=c(models,'Vehicle','NB25'))]
legend_items <- data.table(x=c(1,3.1,5.2,7.1,8.9),category=factor(c(models,'Vehicle','NB25'),levels=c(models,'Vehicle','NB25')))
p6Dleg <- ggplot(legend_items,aes(x,1,fill=category))+geom_point(shape=22,size=5)+
 geom_text(aes(label=ifelse(category=='NB25','NDNB25',as.character(category))),nudge_x=.32,hjust=0,size=10/.pt)+
 scale_fill_manual(values=c(cellcols,treatcols),guide='none')+coord_cartesian(xlim=c(.4,10.5),ylim=c(.7,1.3),clip='off')+
 theme_void()+theme(plot.margin=margin(0,2,0,0))
p6Dheat <- ggplot()+geom_tile(data=hm[structural_NA==FALSE],aes(x,y,fill=z),width=.98,height=.98,color='white',linewidth=.15)+
 geom_tile(data=hm[structural_NA==TRUE],aes(x,y),fill='white',width=.98,height=.98,color='#E5E5E5',linewidth=.15)+
 geom_text(data=hm[structural_NA==TRUE],aes(x,y),label='/',size=11/.pt,color='#404040')+
 scale_fill_gradient2(low=blue,mid='white',high=orange,midpoint=0,limits=c(-2,2),oob=scales::squish,na.value='#D8DADD',name='VST z-score')+
 ggnewscale::new_scale_fill()+geom_tile(data=bars,aes(x,y,fill=category),width=.98,height=.74,color='white',linewidth=.15)+
 scale_fill_manual(values=c(cellcols,treatcols),guide='none')+
 geom_vline(xintercept=c(6.5,12.5),color='white',linewidth=1)+
 scale_x_continuous(breaks=1:18,labels=rep(1:3,6),limits=c(.5,18.5),expand=c(0,0))+
 scale_y_continuous(breaks=c(seq_along(lipgenes),length(lipgenes)+1.2,length(lipgenes)+2.05),
 labels=c(rev(lipgenes),'Treatment','Cell line'),limits=c(.5,length(lipgenes)+2.5),expand=c(0,0))+
 labs(x='Biological replicate',y=NULL)+base_theme+
 theme(panel.grid=element_blank(),axis.text.x=element_text(size=10),axis.text.y=element_text(size=11,face='bold'),
 axis.title.x=element_text(size=12),legend.position='right',legend.title=element_text(size=10),legend.text=element_text(size=10),
 plot.margin=margin(0,4,0,4))
p6D <- (p6Dleg/p6Dheat)+plot_layout(heights=c(.35,3.65))+
 plot_annotation(title='D. Lipogenic genes across biological replicates',
 subtitle='VST-transformed expression; z-scores scaled within each cell line and gene.',
 caption=' / = not applicable because of species-specific gene/ortholog representation.',
 theme=theme(plot.title=element_text(size=16,face='bold'),plot.subtitle=element_text(size=11),plot.title.position='plot',plot.caption=element_text(size=9,hjust=0),plot.margin=margin(3,3,3,3)))
savefig(p6D,'Figure_6','Figure_6D_lipogenic_replicates',7.5,4.55)
# Figure 6E: five archived programs; distinct FDR families explicitly noted.
five <- fread(file.path(data_dir,'five_program_comparison_source.csv'))
want <- c('CUSTOM_SEVEN_GENE_LIPOGENESIS','HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_CHOLESTEROL_HOMEOSTASIS',
 'GOBP_REGULATION_OF_FAT_CELL_DIFFERENTIATION','HALLMARK_ADIPOGENESIS')
names5 <- c('SREBF1-associated lipogenic module','Hallmark fatty-acid metabolism','Hallmark cholesterol homeostasis',
 'GO regulation of fat-cell differentiation','Hallmark adipogenesis')
loc <- c(6.2,5.2,4.2,2.5,1.5)
five <- five[pathway %in% want]
stopifnot(nrow(five)==15)
stopifnot(isTRUE(all.equal(
 five[pathway=='CUSTOM_SEVEN_GENE_LIPOGENESIS'][match(models,cell_line)]$FDR,
 mod[match(models,cell_line)]$BH_q_across_three_prespecified_tests,
 tolerance=1e-12)))
five[,label:=names5[match(pathway,want)]];five[,y0:=loc[match(pathway,want)]]
five[,y:=y0+c(OVCAR5=.19,OVCAR8=0,ID8=-.19)[cell_line]]
five[,cell_line:=factor(cell_line,levels=models)]
five[,sig:=factor(ifelse(FDR<.05,'FDR < 0.05','FDR >= 0.05'),levels=c('FDR < 0.05','FDR >= 0.05'))]
fwrite(five[,.(pathway,label,cell_line,NES,nominal_p,FDR,gene_set_FDR05)],file.path(out,'source_data','lipid_adipogenesis_programs.csv'))
p6E <- ggplot(five,aes(NES,y,color=cell_line,shape=sig))+geom_vline(xintercept=0,color='#AEB4B8',linewidth=.4)+
 geom_point(size=3.2,stroke=1.1)+
 scale_color_manual(values=linecols,breaks=models,name='Cell line')+
 scale_shape_manual(values=c('FDR < 0.05'=16,'FDR >= 0.05'=1),name='Gene-set FDR')+
 scale_y_continuous(breaks=loc,labels=names5,limits=c(1,6.7),expand=c(0,0))+
 scale_x_continuous(limits=c(-3,2.15),breaks=-3:2)+
 labs(title='E. Lipid-metabolic and adipogenesis-related pathway responses to NDNB25',
 subtitle='FDR shown within each gene-set testing family; negative NES is enriched among NDNB25-decreased genes.',
 x='Normalized enrichment score (NES)',y=NULL)+base_theme+
 theme(panel.grid.major.y=element_line(color='#E8EAEC',linewidth=.3),panel.grid.major.x=element_line(color='#EEF0F1',linewidth=.3),
 axis.text.y=element_text(size=10.5,face='bold'),axis.text.x=element_text(size=11),axis.title.x=element_text(size=12),
 legend.position='bottom',legend.box='vertical',legend.title=element_text(size=10,face='bold'),legend.text=element_text(size=10),
 plot.title=element_text(size=14,face='bold'),plot.margin=margin(4,4,4,4))
savefig(p6E,'Figure_6','Figure_6E_five_program_comparison',7.5,4.65)
# Figure 6F: archived predefined sterol-synthesis and SREBF2 effects, no per-gene significance selection.
eff <- fread(file.path(data_dir,'selected_gene_effects_and_mapping.csv'))
sterol <- c('HMGCS1','HMGCR','MVK','IDI1','FDPS','FDFT1','SQLE','DHCR7','DHCR24','SREBF2')
eff <- eff[display_gene %in% sterol]
stopifnot(nrow(eff)==30,all(!is.na(eff$gene_id)),all(!is.na(eff$log2FoldChange)))
eff[,cell_line:=factor(cell_line,levels=models)]
eff[,display_gene:=factor(display_gene,levels=rev(sterol))]
eff[,mark:=ifelse(!is.na(padj)&padj<.05,'*','')]
fwrite(eff[,.(cell_line,display_gene,human_gene,mouse_gene,mapping_method,gene_id,baseMean,log2FC=log2FoldChange,gene_FDR=padj,mark)],
 file.path(out,'source_data','sterol_gene_effects.csv'))
p6F <- ggplot(eff,aes(cell_line,display_gene,fill=log2FoldChange))+geom_tile(color='white',linewidth=.5)+
 geom_text(aes(label=mark),size=12/.pt,fontface='bold')+
 scale_fill_gradient2(low=blue,mid='white',high=orange,midpoint=0,limits=c(-1.5,1.5),oob=scales::squish,na.value='#D8DADD',name='log2FC')+
 scale_x_discrete(position='top')+labs(title='F. Sterol-synthesis and regulatory gene effects',
 subtitle='Paired DESeq2 log2FC (NDNB25 versus vehicle); * gene-level FDR < 0.05.',x=NULL,y=NULL)+
 base_theme+theme(panel.grid=element_blank(),axis.text.x=element_text(size=14,face='bold'),axis.text.y=element_text(size=11,face='bold'),
 legend.position='right',legend.title=element_text(size=10),legend.text=element_text(size=10),plot.margin=margin(4,4,4,4))
savefig(p6F,'Figure_6','Figure_6F_sterol_gene_effects',7.5,4.4)
fig6_linear <- wrap_elements(full=p6AC)/wrap_elements(full=p6D)/wrap_elements(full=p6E)/wrap_elements(full=p6F)+
 plot_layout(heights=c(3.5,4.55,4.65,4.4))
# Review-only export omitted.
p6E_compact <- p6E+
 scale_y_continuous(breaks=loc,labels=c('SREBF1 lipogenic module','Hallmark fatty-acid metabolism','Hallmark cholesterol homeostasis',
 'GO regulation of fat-cell\ndifferentiation','Hallmark adipogenesis'),limits=c(1,6.7),expand=c(0,0))+
 labs(title='E. Lipid and adipogenesis programs',subtitle='NES; filled = FDR < 0.05; open = FDR >= 0.05.')+
 theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=10),
 axis.text.y=element_text(size=10),legend.position='bottom',legend.box='vertical',
 legend.title=element_text(size=9),legend.text=element_text(size=9),plot.margin=margin(4,3,3,3))
p6F_compact <- p6F+
 labs(title='F. Sterol gene effects',subtitle='NDNB25/vehicle; * FDR < 0.05.')+
 guides(fill=guide_colorbar(direction='horizontal',title.position='top',barwidth=grid::unit(1.7,'in'),barheight=grid::unit(.12,'in')))+
 theme(plot.title=element_text(size=14,face='bold'),plot.subtitle=element_text(size=10),
 axis.text.x=element_text(size=9.5,face='bold'),axis.text.y=element_text(size=10),
 legend.position='bottom',legend.title=element_text(size=9),legend.text=element_text(size=9),plot.margin=margin(4,3,3,3))
bottom6 <- (p6E_compact | p6F_compact)+plot_layout(widths=c(3.7,3.8))
fig6 <- wrap_elements(full=p6AC)/wrap_elements(full=p6D)/wrap_elements(full=bottom6)+
 plot_layout(heights=c(2.85,3.5,4.6))
savefig(fig6,'Figure_6','Figure_6_final',7.5,10.95)

