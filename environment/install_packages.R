cran<-c('data.table','ggplot2','patchwork','ggrepel','ggnewscale','ragg','svglite','scales','BiocManager')
missing<-cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
missing<-c('DESeq2','fgsea')[!vapply(c('DESeq2','fgsea'),requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))BiocManager::install(missing,ask=FALSE,update=FALSE)
message('Compare installed versions with package_versions.txt; different versions may change numerical or rendering results.')
