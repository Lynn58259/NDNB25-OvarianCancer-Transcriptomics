if (getRversion() < '4.5.0' || getRversion() >= '4.6.0') stop('Use R 4.5.x (validated: 4.5.1) with Bioconductor 3.22.')
cran<-c('data.table','ggplot2','patchwork','ggrepel','ggnewscale','ragg','svglite','scales','BiocManager')
missing<-cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos='https://cloud.r-project.org')
missing<-c('DESeq2','fgsea')[!vapply(c('DESeq2','fgsea'),requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))BiocManager::install(missing,version='3.22',ask=FALSE,update=FALSE)
message('Compare installed versions with package_versions.txt; different versions may change numerical or rendering results.')
