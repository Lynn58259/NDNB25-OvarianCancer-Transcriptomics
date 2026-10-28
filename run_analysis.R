args<-commandArgs(FALSE);f<-grep('^--file=',args,value=TRUE)
if(length(f))setwd(dirname(normalizePath(sub('^--file=','',f[1]))))
Sys.setlocale('LC_COLLATE','C')
for(d in c('results/quality_control','results/plot_data','figures'))dir.create(d,recursive=TRUE,showWarnings=FALSE)
for(z in c('OVCAR5','OVCAR8','ID8'))for(d in c('differential_expression','normalized_counts','pathway_enrichment'))dir.create(file.path('results',d,z),recursive=TRUE,showWarnings=FALSE)
for(s in sort(list.files('scripts',pattern='^0[1-6]_.*[.]R$',full.names=TRUE))){message('RUN ',s);source(s)}
writeLines(capture.output(sessionInfo()),'environment/sessionInfo.txt')
p<-installed.packages();write.table(p[,c('Package','Version')],'environment/package_versions.txt',sep='\t',quote=FALSE,row.names=FALSE)