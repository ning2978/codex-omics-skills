packages <- c('Seurat','SeuratObject','Matrix','ggplot2','patchwork','jsonlite','data.table','readr','optparse','statmod','edgeR','png','R.utils')
keys <- c('r-seurat','r-seuratobject','r-matrix','r-ggplot2','r-patchwork','r-jsonlite','r-data.table','r-readr','r-optparse','r-statmod','bioconductor-edger','r-png','r-r.utils')
versions <- list('r-base'=as.character(getRversion()))
for (i in seq_along(packages)) {
  suppressPackageStartupMessages(library(packages[i], character.only=TRUE))
  versions[[keys[i]]] <- as.character(packageVersion(packages[i]))
}
cat(jsonlite::toJSON(list(versions=versions, machine=R.version$arch), auto_unbox=TRUE))
