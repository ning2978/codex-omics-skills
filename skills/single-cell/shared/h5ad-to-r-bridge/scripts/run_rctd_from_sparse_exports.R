#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(spacexr)
})

parse_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  args <- list(
    reference_dir = NULL,
    spatial_dir = NULL,
    output_dir = NULL,
    label_key = NULL,
    max_cores = "4",
    test_mode = "FALSE",
    doublet_mode = "doublet",
    n_max_cells = "10000",
    min_umi = "100"
  )
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
  for (required in c("reference_dir", "spatial_dir", "output_dir", "label_key")) {
    if (is.null(args[[required]])) stop("--", gsub("_", "-", required), " is required")
  }
  args$max_cores <- as.integer(args$max_cores)
  args$test_mode <- parse_bool(args$test_mode)
  args$n_max_cells <- as.integer(args$n_max_cells)
  args$min_umi <- as.numeric(args$min_umi)
  args
}

read_sparse_export <- function(input_dir) {
  counts <- as(readMM(gzfile(file.path(input_dir, "counts.mtx.gz"))), "CsparseMatrix")
  features <- fread(file.path(input_dir, "features.tsv.gz"), header = FALSE)
  barcodes <- fread(file.path(input_dir, "barcodes.tsv.gz"), header = FALSE)[[1]]
  rownames(counts) <- make.unique(as.character(features[[2]]))
  colnames(counts) <- as.character(barcodes)
  metadata <- fread(file.path(input_dir, "metadata.csv"), data.table = FALSE)
  rownames(metadata) <- metadata$barcode
  metadata <- metadata[colnames(counts), , drop = FALSE]
  if (!identical(rownames(metadata), colnames(counts))) stop("metadata barcode order mismatch for ", input_dir)
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

make_reference <- function(reference_export, args) {
  if (!args$label_key %in% colnames(reference_export$metadata)) {
    stop("Reference metadata missing label key: ", args$label_key)
  }
  cell_types <- sanitize_cell_types(reference_export$metadata[[args$label_key]], args$output_dir)
  Reference(
    reference_export$counts,
    factor(cell_types),
    nUMI = Matrix::colSums(reference_export$counts),
    require_int = TRUE,
    n_max_cells = args$n_max_cells,
    min_UMI = args$min_umi
  )
}

make_spatial <- function(spatial_export, spatial_dir) {
  coords <- fread(file.path(spatial_dir, "spatial", "tissue_positions.csv"), data.table = FALSE)
  rownames(coords) <- coords$barcode
  coords <- coords[colnames(spatial_export$counts), c("x", "y"), drop = FALSE]
  if (!identical(rownames(coords), colnames(spatial_export$counts))) stop("spatial coordinate order mismatch")
  SpatialRNA(coords, spatial_export$counts, nUMI = Matrix::colSums(spatial_export$counts), require_int = TRUE)
}

extract_weights <- function(rctd_obj) {
  if (!is.null(rctd_obj@results$weights)) return(as.matrix(rctd_obj@results$weights))
  if (!is.null(rctd_obj@results$weights_doublet)) return(as.matrix(rctd_obj@results$weights_doublet))
  NULL
}

main <- function() {
  args <- parse_args()
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  reference_export <- read_sparse_export(args$reference_dir)
  spatial_export <- read_sparse_export(args$spatial_dir)

  common_genes <- intersect(rownames(reference_export$counts), rownames(spatial_export$counts))
  if (length(common_genes) == 0) stop("No shared genes between reference and spatial exports")
  reference_export$counts <- reference_export$counts[common_genes, , drop = FALSE]
  spatial_export$counts <- spatial_export$counts[common_genes, , drop = FALSE]

  reference <- make_reference(reference_export, args)
  spatial <- make_spatial(spatial_export, args$spatial_dir)
  saveRDS(reference, file.path(args$output_dir, "reference.rds"))
  saveRDS(spatial, file.path(args$output_dir, "spatial_rna.rds"))

  rctd <- create.RCTD(spatial, reference, max_cores = args$max_cores, test_mode = args$test_mode)
  rctd <- run.RCTD(rctd, doublet_mode = args$doublet_mode)
  saveRDS(rctd, file.path(args$output_dir, "rctd_result.rds"))

  weights <- extract_weights(rctd)
  if (!is.null(weights)) {
    weights_df <- data.frame(barcode = rownames(weights), weights, check.names = FALSE)
    fwrite(weights_df, file.path(args$output_dir, "rctd_weights.csv"))
  }

  manifest <- data.frame(
    metric = c(
      "reference_dir", "spatial_dir", "shared_genes", "reference_cells",
      "spatial_locations", "max_cores", "test_mode", "doublet_mode"
    ),
    value = c(
      args$reference_dir, args$spatial_dir, length(common_genes), ncol(reference_export$counts),
      ncol(spatial_export$counts), args$max_cores, args$test_mode, args$doublet_mode
    )
  )
  write.csv(manifest, file.path(args$output_dir, "run_manifest.csv"), row.names = FALSE)
  cat("Wrote RCTD result directory:", args$output_dir, "\n")
}

main()
