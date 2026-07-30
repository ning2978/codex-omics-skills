#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(spacexr)
})

parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  args <- list(input_dir = NULL, type = NULL, output_rds = NULL, label_key = NULL, n_max_cells = "10000", min_umi = "100")
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
  if (is.null(args$type) || !args$type %in% c("reference", "spatial")) stop("--type must be reference or spatial")
  if (is.null(args$output_rds)) stop("--output-rds is required")
  args$n_max_cells <- as.integer(args$n_max_cells)
  args$min_umi <- as.numeric(args$min_umi)
  args
}

read_sparse_export <- function(input_dir) {
  counts <- as(readMM(gzfile(file.path(input_dir, "counts.mtx.gz"))), "CsparseMatrix")
  features <- fread(file.path(input_dir, "features.tsv.gz"), header = FALSE)
  barcodes <- fread(file.path(input_dir, "barcodes.tsv.gz"), header = FALSE)[[1]]
  if (nrow(features) != nrow(counts)) stop("features rows do not match counts rows")
  if (length(barcodes) != ncol(counts)) stop("barcodes do not match counts columns")
  rownames(counts) <- make.unique(as.character(features[[2]]))
  colnames(counts) <- as.character(barcodes)
  metadata <- fread(file.path(input_dir, "metadata.csv"), data.table = FALSE)
  if (!"barcode" %in% colnames(metadata)) stop("metadata.csv must contain barcode")
  rownames(metadata) <- metadata$barcode
  metadata <- metadata[colnames(counts), , drop = FALSE]
  if (!identical(rownames(metadata), colnames(counts))) stop("metadata barcode order mismatch")
  list(counts = counts, metadata = metadata)
}

sanitize_cell_types <- function(cell_types, output_dir) {
  original <- as.character(cell_types)
  sanitized <- gsub("[^A-Za-z0-9_.-]+", "_", original)
  sanitized <- gsub("_+", "_", sanitized)
  sanitized <- gsub("^_|_$", "", sanitized)
  mapping <- unique(data.frame(original_cell_type = original, rctd_cell_type = sanitized))
  write.csv(mapping, file.path(output_dir, "cell_type_label_mapping.csv"), row.names = FALSE)
  sanitized
}

build_reference <- function(export, args, output_dir) {
  if (is.null(args$label_key) || !args$label_key %in% colnames(export$metadata)) {
    stop("--label-key must name a metadata column for reference exports")
  }
  cell_types <- sanitize_cell_types(export$metadata[[args$label_key]], output_dir)
  nUMI <- Matrix::colSums(export$counts)
  Reference(
    export$counts,
    factor(cell_types),
    nUMI = nUMI,
    require_int = TRUE,
    n_max_cells = args$n_max_cells,
    min_UMI = args$min_umi
  )
}

build_spatial <- function(export, input_dir) {
  spatial_path <- file.path(input_dir, "spatial", "tissue_positions.csv")
  if (!file.exists(spatial_path)) stop("Missing spatial/tissue_positions.csv for spatial export")
  coords <- fread(spatial_path, data.table = FALSE)
  if (!all(c("barcode", "x", "y") %in% colnames(coords))) stop("tissue_positions.csv must contain barcode, x, y")
  rownames(coords) <- coords$barcode
  coords <- coords[colnames(export$counts), c("x", "y"), drop = FALSE]
  if (!identical(rownames(coords), colnames(export$counts))) stop("Spatial coordinate order mismatch")
  nUMI <- Matrix::colSums(export$counts)
  SpatialRNA(coords, export$counts, nUMI = nUMI, require_int = TRUE)
}

main <- function() {
  args <- parse_args()
  dir.create(dirname(args$output_rds), recursive = TRUE, showWarnings = FALSE)
  export <- read_sparse_export(args$input_dir)
  output_dir <- dirname(args$output_rds)
  object <- if (args$type == "reference") {
    build_reference(export, args, output_dir)
  } else {
    build_spatial(export, args$input_dir)
  }
  saveRDS(object, args$output_rds)
  validation <- data.frame(
    metric = c("type", "n_genes", "n_observations", "counts_nnz"),
    value = c(args$type, nrow(export$counts), ncol(export$counts), length(export$counts@x))
  )
  write.csv(validation, paste0(args$output_rds, ".validation.csv"), row.names = FALSE)
  cat("Wrote RCTD input object:", args$output_rds, "\n")
}

main()
