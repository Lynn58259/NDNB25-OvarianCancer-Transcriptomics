suppressPackageStartupMessages({library(data.table);library(DESeq2)})
models<-c('OVCAR5','OVCAR8','ID8');meta<-fread('data/metadata/samples.csv',data.table=FALSE)
raw<-read.delim('data/counts/provider_expression_matrix.tsv',check.names=FALSE)
stopifnot(nrow(meta)==18,!anyDuplicated(meta$sample_id),!anyDuplicated(raw$gene_id))
inputs<-list();pre<-list()
for(z in models){
 mm<-meta[meta$cell_line==z,];rownames(mm)<-mm$sample_id
 stopifnot(nrow(mm)==6,all(table(mm$treatment)==3),all(table(mm$batch,mm$treatment)==1))
 mm$treatment<-factor(mm$treatment,levels=c('CON','NB25'));mm$batch<-factor(mm$batch)
 prefix<-if(z=='ID8')'^ENSMUSG' else '^ENSG';stopifnot(all(mm$species==if(z=='ID8')'mouse' else 'human'))
 aa<-raw[grepl(prefix,raw$gene_id),c('gene_id','gene_name','gene_biotype')]
 cc<-as.matrix(raw[grepl(prefix,raw$gene_id),paste0(mm$sample_id,'_count')]);rownames(cc)<-aa$gene_id;colnames(cc)<-mm$sample_id
 stopifnot(all(is.finite(cc)),all(cc>=0));rounded<-round(cc);storage.mode(rounded)<-'integer';keep<-rowSums(rounded>=10)>=3
 inputs[[z]]<-list(metadata=mm,annotation=aa[keep,],counts=rounded[keep,,drop=FALSE],original=cc[keep,,drop=FALSE])
 pre[[z]]<-data.frame(cell_line=z,samples=6,genes_before=nrow(cc),genes_after=sum(keep),fractional_entries=sum(abs(cc-rounded)>1e-7))
}
saveRDS(inputs,'results/quality_control/preprocessed_inputs.rds');fwrite(rbindlist(pre),'results/quality_control/preprocessing_summary.csv')
