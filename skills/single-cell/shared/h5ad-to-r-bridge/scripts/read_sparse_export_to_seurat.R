#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
})

parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  args <- list(input_dir = NULL, output_rds = NULL, assay = "RNA", project = "h5ad_bridge")
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--") || i == length(argv)) stop("Invalid argument near: ", key)
    value <- argv[[i + 1]]
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    args[[name]] <- value
    i <- i + 2
  }
  if (is.null(args$input_dir)) stop("--input-dir is required")
  if (is.null(args$output_rds)) stop("--output-rds is required")
  args
}

read_vector_gz <- function(path) {
  fread(path, header = FALSE)[[1]]
}

read_sparse_export <- function(input_dir) {
  counts <- as(readMM(gzfile(file.path(input_dir, "counts.mtx.gz"))), "CsparseMatrix")
  features <- fread(file.path(input_dir, "features.tsv.gz"), header = FALSE)
  barcodes <- read_vector_gz(file.path(input_dir, "barcodes.tsv.gz"))
  if (nrow(features) != nrow(counts)) stop("features rows do not match counts rows")
  if (length(barcodes) != ncol(counts)) stop("barcodes do not match counts columns")
  gene_names <- make.unique(as.character(features[[2]]))
  rownames(counts) <- gene_names
  colnames(counts) <- as.character(barcodes)
  metadata <- fread(file.path(input_dir, "metadata.csv"), data.table = FALSE)
  if (!"barcode" %in% colnames(metadata)) stop("metadata.csv must contain barcode")
  rownames(metadata) <- metadata$barcode
  metadata <- metadata[colnames(counts), , drop = FALSE]
  if (!identical(rownames(metadata), colnames(counts))) stop("metadata barcode order mismatch")
  list(counts = counts, metadata = metadata)
}

main <- function() {
  args <- parse_args()
  export <- read_sparse_export(args$input_dir)
  obj <- CreateSeuratObject(
    counts = export$counts,
    assay = args$assay,
    project = args$project,
    meta.data = export$metadata
  )
  spatial_path <- file.path(args$input_dir, "spatial", "tissue_positions.csv")
  if (file.exists(spatial_path)) {
    coords <- fread(spatial_path, data.table = FALSE)
    rownames(coords) <- coords$barcode
    coords <- coords[colnames(obj), , drop = FALSE]
    obj[["h5ad_bridge_x"]] <- coords$x
    obj[["h5ad_bridge_y"]] <- coords$y
    obj@misc$h5ad_bridge_spatial_coordinates <- coords
  }
  dir.create(dirname(args$output_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, args$output_rds)
  validation <- data.frame(
    metric = c("n_features", "n_observations", "metadata_rows", "has_spatial_coords"),
    value = c(nrow(obj), ncol(obj), nrow(obj@meta.data), file.exists(spatial_path))
  )
  write.csv(validation, paste0(args$output_rds, ".validation.csv"), row.names = FALSE)
  cat("Wrote Seurat object:", args$output_rds, "\n")
}

main()
