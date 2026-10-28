suppressPackageStartupMessages(library(fgsea))
get_stats<-function(z){d<-fread(file.path('results/plot_data',paste0(z,'_ranked_GSEA.csv')));d<-d[!is.na(gene_label)&gene_label!=''&!is.na(stat)];d[,key:=toupper(gene_label)];d<-d[!duplicated(key)];st<-d$stat;names(st)<-d$key;sort(st,decreasing=TRUE)}
for(i in seq_along(models)){
 z<-models[i];st<-get_stats(z)
 for(family in c('Hallmark','GOBP')){
  message('GSEA: ',z,' ',family);sets<-readRDS(file.path('data/gene_sets',paste0(if(z=='ID8')'mouse' else 'human','_',family,'.rds')))
  set.seed(20260922+i);f<-as.data.table(fgseaMultilevel(pathways=sets,stats=st,minSize=10,maxSize=500,eps=0,nproc=1))
  f[,leadingEdge:=vapply(leadingEdge,paste,collapse=';',FUN.VALUE=character(1))];f[,cell_line:=z]
  fwrite(f,file.path('results/pathway_enrichment',z,paste0(family,'_GSEA.csv')))
  if(family=='Hallmark')fwrite(f[order(padj,-abs(NES))],file.path('results/plot_data',paste0(z,'_Hallmark_only_seeded_multilevel_GSEA.csv')))
 }
}
