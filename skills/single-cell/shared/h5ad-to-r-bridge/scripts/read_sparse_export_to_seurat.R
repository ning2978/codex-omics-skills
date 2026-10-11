#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
  library(jsonlite)
})
parse_args <- function() {
  argv <- commandArgs(trailingOnly=TRUE)
  args <- list(input_dir=NULL, output_rds=NULL, assay="RNA", project="h5ad_bridge", overwrite="FALSE")
  if (length(argv) %% 2 != 0) stop("Arguments require --key value pairs")
  if (length(argv)) for(i in seq(1, length(argv), by=2)) {
    name <- gsub("-", "_", sub("^--", "", argv[[i]]))
    if (!name %in% names(args)) stop("Unknown argument: ", argv[[i]])
    args[[name]] <- argv[[i+1]]
  }
  if (is.null(args$input_dir) || is.null(args$output_rds)) stop("--input-dir and --output-rds required")
  args
}
read_keyed <- function(path, keys) {
  tab <- fread(path, data.table=FALSE, colClasses=c(barcode="character"))
  if (!"barcode" %in% names(tab) || anyNA(tab$barcode) || anyDuplicated(tab$barcode)) stop("Invalid/duplicate barcode in ", path)
  if (!setequal(tab$barcode, keys)) stop("Cell key set mismatch in ", path)
  rownames(tab) <- tab$barcode
  tab[keys, , drop=FALSE]
}
main <- function() {
  args <- parse_args()
  if (file.exists(args$output_rds) && !isTRUE(as.logical(args$overwrite))) stop("Output exists; choose another path or --overwrite TRUE")
  counts <- as(readMM(gzfile(file.path(args$input_dir,"counts.mtx.gz"))), "CsparseMatrix")
  features <- fread(file.path(args$input_dir,"features.tsv.gz"), header=FALSE, colClasses="character", data.table=FALSE)
  barcodes <- fread(file.path(args$input_dir,"barcodes.tsv.gz"), header=FALSE, colClasses="character")[[1]]
  if (nrow(features)!=nrow(counts) || length(barcodes)!=ncol(counts) || ncol(features)<2) stop("Feature/barcode dimensions mismatch")
  if (anyDuplicated(barcodes) || anyNA(barcodes) || any(barcodes=="")) stop("Invalid cell keys")
  if (anyDuplicated(features[[1]]) || anyNA(features[[1]]) || any(features[[1]]=="")) stop("Unique feature IDs required")
  if (any(!is.finite(counts@x)) || any(counts@x<0) || any(abs(counts@x-round(counts@x))>1e-6)) stop("Invalid original raw counts")
  # Seurat disallows underscores/pipes in feature names. Preserve reversible IDs.
  seurat_features <- make.unique(gsub("[_|]", "-", features[[1]]))
  rownames(counts) <- seurat_features
  colnames(counts) <- barcodes
  meta <- read_keyed(file.path(args$input_dir,"metadata.csv"), barcodes)
  fmap <- data.frame(source_feature_id=features[[1]], gene_symbol=features[[2]],
                     seurat_feature=seurat_features, source_row=seq_len(nrow(features)), stringsAsFactors=FALSE)
  feature_path <- file.path(args$input_dir,"feature_metadata.csv")
  if (file.exists(feature_path)) {
    original <- fread(feature_path, data.table=FALSE)
    if (nrow(original)!=nrow(fmap)) stop("Feature metadata dimension mismatch")
    fmap <- cbind(fmap, original[,setdiff(names(original),names(fmap)),drop=FALSE])
  }
  obj <- CreateSeuratObject(counts=counts, assay=args$assay, project=args$project, meta.data=meta)
  if (!identical(colnames(obj),barcodes) || !identical(rownames(obj[[args$assay]]),seurat_features)) stop("Seurat altered keys/order")
  rownames(fmap) <- seurat_features
  obj[[args$assay]] <- AddMetaData(obj[[args$assay]], fmap)
  obj@misc$h5ad_bridge_feature_mapping <- fmap
  manifest_path <- file.path(args$input_dir,"source_manifest.json")
  if (!file.exists(manifest_path)) stop("Source manifest required to declare reductions")
  manifest <- fromJSON(manifest_path, simplifyVector=TRUE)
  obj@misc$h5ad_bridge_source_manifest <- manifest
  reduction_map <- data.frame(source_file=character(), seurat_reduction=character(), dimensions=integer())
  paths <- unlist(manifest$reductions_exported, use.names=FALSE)
  used <- character()
  for (i in seq_along(paths)) {
    relative <- paths[[i]]
    if (!grepl("^reductions/[A-Za-z0-9_-]+\\.csv$", relative)) stop("Invalid reduction manifest path")
    tab <- read_keyed(file.path(args$input_dir,relative),barcodes)
    values <- as.matrix(tab[,setdiff(names(tab),"barcode"),drop=FALSE])
    if (!is.numeric(values) || !all(is.finite(values)) || ncol(values)<1) stop("Invalid reduction values")
    stem <- sub("\\.csv$", "", basename(relative))
    name <- make.unique(c(used,tolower(sub("^X_","",stem))))[length(used)+1]
    used <- c(used,name)
    original_columns <- colnames(values)
    reduction_key <- paste0("BR",i,"_")
    colnames(values) <- paste0(reduction_key,seq_len(ncol(values)))
    obj[[name]] <- CreateDimReducObject(embeddings=values, key=reduction_key, assay=args$assay)
    if (!isTRUE(all.equal(unname(Embeddings(obj[[name]])),unname(values),tolerance=1e-12))) stop("Reduction numeric mismatch")
    obj@misc[[paste0("h5ad_bridge_reduction_",name)]] <- list(source_file=relative, source_columns=original_columns)
    reduction_map <- rbind(reduction_map,data.frame(source_file=relative,seurat_reduction=name,dimensions=ncol(values)))
  }
  spatial_path <- file.path(args$input_dir,"spatial","tissue_positions.csv")
  if (file.exists(spatial_path)) {
    coords <- read_keyed(spatial_path,barcodes)
    if (!all(c("x","y") %in% names(coords)) || !all(is.finite(as.matrix(coords[,c("x","y")])))) stop("Invalid spatial coordinates")
    obj[["h5ad_bridge_x"]] <- coords$x
    obj[["h5ad_bridge_y"]] <- coords$y
    obj@misc$h5ad_bridge_spatial_coordinates <- coords
  }
  observed <- GetAssayData(obj,assay=args$assay,layer="counts")
  if (!isTRUE(all.equal(as(observed,"CsparseMatrix"),counts,tolerance=0))) stop("Seurat count parity failed")
  obj@misc$h5ad_bridge_reduction_mapping <- reduction_map
  obj@misc$h5ad_bridge_runtime <- list(R=R.version.string,libPaths=.libPaths(),
                                     Seurat=as.character(packageVersion("Seurat")),package_path=find.package("Seurat"))
  dir.create(dirname(args$output_rds),recursive=TRUE,showWarnings=FALSE)
  saveRDS(obj,args$output_rds)
  reopened <- readRDS(args$output_rds)
  if (!identical(colnames(reopened),barcodes) || length(Reductions(reopened))!=length(paths)) stop("RDS reopen parity failed")
  write.csv(fmap,paste0(args$output_rds,".feature_mapping.csv"),row.names=FALSE)
  write.csv(reduction_map,paste0(args$output_rds,".reduction_mapping.csv"),row.names=FALSE)
  validation <- data.frame(metric=c("n_features","n_observations","count_parity","cell_order","feature_mapping","reductions_preserved","has_spatial_coords"),
                           value=c(nrow(obj),ncol(obj),TRUE,TRUE,TRUE,nrow(reduction_map),file.exists(spatial_path)))
  write.csv(validation,paste0(args$output_rds,".validation.csv"),row.names=FALSE)
  cat("Wrote Seurat object:",args$output_rds,"\n")
}
main()
