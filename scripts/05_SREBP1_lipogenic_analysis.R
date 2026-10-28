module_sets<-readRDS('data/gene_sets/SREBF1_module.rds');module<-list();fx<-list();hm<-list();five<-list();selection<-fread('data/gene_sets/selected_gene_mapping.csv')
for(i in seq_along(models)){
 z<-models[i];st<-get_stats(z);gs<-module_sets[[z]];set.seed(20260922+i)
 f<-fgseaMultilevel(list(PRESPECIFIED_SREBF1_ASSOCIATED_LIPOGENIC_MODULE=gs),st,minSize=5,maxSize=30,eps=0,nproc=1)
 module[[z]]<-data.table(cell_line=z,gene_set_size=f$size,ES=f$ES,NES=f$NES,pvalue=f$pval,log2err=f$log2err,leading_edge=paste(f$leadingEdge[[1]],collapse=';'),present_genes=paste(intersect(gs,names(st)),collapse=';'),mapping_count=length(gs),tested_gene_count=length(st))
 de<-fread(file.path('results/differential_expression',z,'all_genes.csv'));sel<-selection[cell_line==z];sel[,gene_name:=NULL];q<-merge(sel,de,by='gene_id',all.x=TRUE,sort=FALSE)
 structural<-q$display_gene=='SCD (Scd2)' & z!='ID8'
 stopifnot(all(!is.na(q$baseMean[!structural])),all(grepl(if(z=='ID8')'^ENSMUSG' else '^ENSG',q$gene_id[!structural])));fx[[z]]<-q
 v<-fread(file.path('results/differential_expression',z,'VST.csv'));vv<-melt(v,id.vars='gene_id',variable.name='sample_id',value.name='vst');vv<-merge(vv,q[,.(gene_id,display_gene,arm)],by='gene_id',allow.cartesian=TRUE)
 vv[,cell_line:=z];vv[,z:=as.numeric(scale(vst)),by=.(display_gene)];hm[[z]]<-vv
 h<-fread(file.path('results/pathway_enrichment',z,'Hallmark_GSEA.csv'))[pathway%in%c('HALLMARK_FATTY_ACID_METABOLISM','HALLMARK_CHOLESTEROL_HOMEOSTASIS','HALLMARK_ADIPOGENESIS')]
 g<-fread(file.path('results/pathway_enrichment',z,'GOBP_GSEA.csv'))[pathway=='GOBP_REGULATION_OF_FAT_CELL_DIFFERENTIATION']
 five[[z]]<-rbindlist(list(h,g),fill=TRUE)[,.(pathway,cell_line,NES,nominal_p=pval,FDR=padj,gene_set_size=size)]
}
mod<-rbindlist(module);mod[,BH_q_across_three_prespecified_tests:=p.adjust(pvalue,'BH')];fwrite(mod,'results/plot_data/SREBP1_lipogenic_module_statistics.csv')
fwrite(rbindlist(fx,fill=TRUE),'results/plot_data/selected_gene_effects_and_mapping.csv');fwrite(rbindlist(hm,fill=TRUE),'results/plot_data/selected_gene_VST_by_sample.csv')
five<-rbindlist(c(five,list(mod[,.(pathway='CUSTOM_SEVEN_GENE_LIPOGENESIS',cell_line,NES,nominal_p=pvalue,FDR=BH_q_across_three_prespecified_tests,gene_set_size)])),fill=TRUE)
setorder(five,pathway,cell_line);five[,gene_set_FDR05:=FDR<.05];fwrite(five,'results/plot_data/five_program_comparison_source.csv');fwrite(as.data.table(meta),'results/plot_data/sample_metadata_confirmed_pairs.csv')

for(model_name in models)fwrite(mod[cell_line==model_name],file.path('results/pathway_enrichment',model_name,'SREBP1_lipogenic_module_GSEA.csv'))
