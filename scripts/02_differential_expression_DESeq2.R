inputs<-readRDS('results/quality_control/preprocessed_inputs.rds');summaries<-list()
for(z in models){
 message('DESeq2: ',z);x<-inputs[[z]];dds<-DESeqDataSetFromMatrix(x$counts,x$metadata,~batch+treatment)
 dds<-DESeq(dds,quiet=TRUE);res<-results(dds,contrast=c('treatment','NB25','CON'),alpha=.05)
 rr<-cbind(x$annotation,as.data.frame(res));rownames(rr)<-NULL
 rr$fractional_entries_original<-rowSums(abs(x$original-x$counts)>1e-7);rr$max_cooks_distance<-apply(assays(dds)[['cooks']],1,max,na.rm=TRUE);rr$converged<-mcols(dds)$betaConv
 rr$FDR05<-!is.na(rr$padj)&rr$padj<.05;rr$FDR05_absLFC1<-rr$FDR05&abs(rr$log2FoldChange)>=1;rr<-rr[order(rr$padj,rr$pvalue,na.last=TRUE),]
 dest<-file.path('results/differential_expression',z);write.csv(rr,file.path(dest,'all_genes.csv'),row.names=FALSE);write.csv(rr[rr$FDR05,],file.path(dest,'DEGs_FDR05.csv'),row.names=FALSE)
 norm<-counts(dds,normalized=TRUE);v<-assay(varianceStabilizingTransformation(dds,blind=FALSE))
 write.csv(data.frame(gene_id=rownames(norm),norm,check.names=FALSE),file.path('results/normalized_counts',z,'normalized_counts.csv'),row.names=FALSE)
 write.csv(data.frame(gene_id=rownames(v),v,check.names=FALSE),file.path(dest,'VST.csv'),row.names=FALSE);saveRDS(dds,file.path(dest,'dds.rds'))
 d<-as.data.table(rr)[!is.na(padj)];d[,gene_label:=fifelse(!is.na(gene_name)&gene_name!='',gene_name,gene_id)];setorder(d,padj,pvalue,gene_id);setorder(d,-stat)
 fwrite(d[,.(gene_id,gene_label,stat,log2FoldChange,pvalue,padj)],file.path('results/plot_data',paste0(z,'_ranked_GSEA.csv')))
 summaries[[z]]<-data.frame(cell_line=z,design='~ batch + treatment',reference='CON (Vehicle)',contrast='NDNB25 versus Vehicle',DEGs=sum(rr$FDR05),nonconverged=sum(!rr$converged))
}
fwrite(rbindlist(summaries),'results/quality_control/DE_summary.csv')
