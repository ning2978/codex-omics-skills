#!/usr/bin/env Rscript

# Render one native Seurat DotPlot per declared grouping.  The unmodified
# plot$data table is written before the fixed display-only colour scale and
# patchwork marker-program strip are added.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10L) {
  stop(
    paste(
      "usage: seurat_dotplot_renderer.R counts.mtx data.mtx marker_contract.csv",
      "cells.csv group_levels.csv native_dir figure_dir width height dpi"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))

counts_path <- args[[1L]]
data_path <- args[[2L]]
marker_contract_path <- args[[3L]]
cells_path <- args[[4L]]
group_levels_path <- args[[5L]]
native_dir <- args[[6L]]
figure_dir <- args[[7L]]
figure_width <- as.numeric(args[[8L]])
figure_height <- as.numeric(args[[9L]])
figure_dpi <- as.numeric(args[[10L]])

if (
  !is.finite(figure_width) || figure_width <= 0 ||
  !is.finite(figure_height) || figure_height <= 0 ||
  !is.finite(figure_dpi) || figure_dpi <= 0
) {
  stop("figure width, height, and dpi must be positive", call. = FALSE)
}

marker_contract <- read.csv(
  marker_contract_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cells <- read.csv(cells_path, stringsAsFactors = FALSE, check.names = FALSE)
group_levels <- read.csv(
  group_levels_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_marker_columns <- c(
  "gene",
  "marker_program",
  "marker_program_short_label",
  "marker_program_color",
  "marker_program_size",
  "marker_role",
  "program_order",
  "program_gene_order",
  "marker_order"
)
required_level_columns <- c(
  "grouping_level",
  "group_order",
  "group_label",
  "group_title"
)
if (!all(required_marker_columns %in% colnames(marker_contract))) {
  stop("marker annotation-strip contract is incomplete", call. = FALSE)
}
if (!all(required_level_columns %in% colnames(group_levels))) {
  stop("group-level contract is incomplete", call. = FALSE)
}
if (!("obs_name" %in% colnames(cells))) {
  stop("cells.csv is missing obs_name", call. = FALSE)
}

markers <- marker_contract$gene
groupings <- unique(group_levels$grouping_level)
if (anyDuplicated(markers) || anyDuplicated(cells$obs_name)) {
  stop("markers and obs_name must be unique", call. = FALSE)
}
if (!all(groupings %in% colnames(cells))) {
  stop("cells.csv is missing one or more grouping columns", call. = FALSE)
}
if (!identical(marker_contract$marker_order, seq_len(nrow(marker_contract)) - 1L)) {
  stop("marker order must be zero-based and contiguous", call. = FALSE)
}

counts <- as(readMM(counts_path), "dgCMatrix")
data <- as(readMM(data_path), "dgCMatrix")
if (!identical(dim(counts), c(length(markers), nrow(cells)))) {
  stop("counts matrix dimensions changed", call. = FALSE)
}
if (!identical(dim(data), dim(counts))) {
  stop("data and counts matrix dimensions differ", call. = FALSE)
}
if (length(counts@x) && (any(!is.finite(counts@x)) || any(counts@x < 0))) {
  stop("counts matrix contains invalid values", call. = FALSE)
}
if (length(data@x) && (any(!is.finite(data@x)) || any(data@x < 0))) {
  stop("data matrix contains invalid values", call. = FALSE)
}

rownames(counts) <- markers
colnames(counts) <- cells$obs_name
rownames(data) <- markers
colnames(data) <- cells$obs_name
meta <- data.frame(row.names = cells$obs_name, check.names = FALSE)
for (grouping in groupings) {
  levels_for_grouping <- group_levels$group_label[
    group_levels$grouping_level == grouping
  ]
  if (!setequal(unique(cells[[grouping]]), levels_for_grouping)) {
    stop(sprintf("%s identities do not match declared levels", grouping), call. = FALSE)
  }
  meta[[grouping]] <- factor(
    cells[[grouping]],
    levels = levels_for_grouping,
    ordered = TRUE
  )
}

object <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  meta.data = meta,
  project = "native_seurat_dotplot_renderer"
)
LayerData(object = object[["RNA"]], layer = "data") <- data
DefaultAssay(object) <- "RNA"

dir.create(native_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

program_rows <- !duplicated(marker_contract$program_order)
strip_data <- marker_contract[program_rows, c(
  "marker_program",
  "marker_program_short_label",
  "marker_program_color",
  "marker_program_size",
  "marker_order"
)]
strip_data$xmin <- strip_data$marker_order + 0.5
strip_data$xmax <- strip_data$xmin + strip_data$marker_program_size
strip_data$x <- (strip_data$xmin + strip_data$xmax) / 2
strip_data$ymin <- 0
strip_data$ymax <- 1
strip_data$y <- 0.5

make_strip <- function(style) ggplot2::ggplot(strip_data) +
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = marker_program_color
    ),
    color = "white",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = x, y = y, label = marker_program_short_label),
    color = "black",
    fontface = "bold",
    size = style$strip_font_pt / ggplot2::.pt,
    lineheight = 0.9
  ) +
  ggplot2::scale_fill_identity(guide = "none") +
  ggplot2::scale_x_continuous(
    limits = c(0.5, length(markers) + 0.5),
    expand = ggplot2::expansion(add = 0)
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

for (grouping in groupings) {
  levels_for_grouping <- group_levels$group_label[
    group_levels$grouping_level == grouping
  ]
  title_values <- unique(group_levels$group_title[
    group_levels$grouping_level == grouping
  ])
  if (length(title_values) != 1L) {
    stop(sprintf("%s has an invalid title contract", grouping), call. = FALSE)
  }
  Idents(object) <- factor(
    cells[[grouping]],
    levels = levels_for_grouping,
    ordered = TRUE
  )

  # This exact native call is the sole source of DotPlot statistics and both
  # native legends.  Export plot$data before any rendering override.
  plot <- DotPlot(object, features=markers, assay="RNA")
  write.csv(
    plot$data,
    file = file.path(native_dir, paste0(grouping, "_seurat_native_pdata.csv")),
    row.names = FALSE,
    quote = TRUE,
    na = "NA"
  )

  style <- as.list(group_levels[group_levels$grouping_level == grouping, ][1, ])
  figure_height <- style$height_inches
  render_plot <- plot +
    ggplot2::scale_x_discrete(limits = markers, expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC",
      mid = "#FFFFFF",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = if (length(levels_for_grouping) == 1) "Expression (unscaled)" else "Average Expression"
    ) +
    ggplot2::labs(x = NULL, y = NULL, colour = if (length(levels_for_grouping) == 1) "Expression (unscaled)" else "Average Expression") +
    ggplot2::guides(colour = ggplot2::guide_colourbar(barheight = grid::unit(0.65, "in"))) +
    ggplot2::theme_classic(base_size = style$group_font_pt) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 60,
        hjust = 1,
        vjust = 1,
        size = style$gene_font_pt
      ),
      axis.text.y = ggplot2::element_text(size = style$group_font_pt, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.3),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      legend.key.height = grid::unit(0.18, "in"),
      legend.spacing.y = grid::unit(2, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.title = ggplot2::element_text(size = style$legend_font_pt),
      legend.text = ggplot2::element_text(size = style$legend_font_pt)
    )

  panels <- make_strip(style) / render_plot +
    patchwork::plot_layout(heights = grid::unit(c(style$strip_height_inches, 1), c("in", "null")))
  combined_plot <- (panels | patchwork::guide_area()) +
    patchwork::plot_layout(widths = grid::unit(c(1, 2.2 * style$legend_font_pt / 11), c("null", "in")), guides = "collect") +
    patchwork::plot_annotation(
      title = title_values[[1L]],
      subtitle = if (is.na(style$group_subtitle)) NULL else style$group_subtitle,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = style$title_font_pt, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = style$group_font_pt)))

  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0(grouping, "_seurat_dotplot.png")),
    plot = combined_plot,
    width = figure_width,
    height = figure_height,
    units = "in",
    dpi = figure_dpi,
    bg = "white"
  )
}

version_record <- data.frame(
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  seurat_version = as.character(packageVersion("Seurat")),
  seuratobject_version = as.character(packageVersion("SeuratObject")),
  patchwork_version = as.character(packageVersion("patchwork")),
  dotplot_call = 'DotPlot(object, features=markers, assay="RNA")',
  render_palette = "#2166AC->#FFFFFF->#B2182B",
  render_midpoint = 0,
  render_limit_min = -2.5,
  render_limit_max = 2.5,
  annotation_strip_legend = FALSE,
  x_axis_text_angle_degrees = 60,
  x_axis_text_hjust = 1,
  x_axis_text_vjust = 1,
  x_axis_text_font_size_pt = style$gene_font_pt,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(
  version_record,
  file = file.path(native_dir, "seurat_render_environment.csv"),
  row.names = FALSE,
  quote = TRUE
)

cat(
  sprintf(
    paste(
      "PASS native Seurat DotPlot renderer: R=%s Seurat=%s",
      "SeuratObject=%s groupings=%d markers=%d\n"
    ),
    version_record$r_version[[1L]],
    version_record$seurat_version[[1L]],
    version_record$seuratobject_version[[1L]],
    length(groupings),
    length(markers)
  )
)
