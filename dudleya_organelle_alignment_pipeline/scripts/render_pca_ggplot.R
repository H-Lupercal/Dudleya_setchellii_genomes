#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "usage: render_pca_ggplot.R COORDINATES_TSV VARIANCE_TSV OUTPUT_PREFIX ORGANELLE",
    call. = FALSE
  )
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
source(file.path(dirname(normalizePath(script_arg)), "dudleya_plotting.R"))

coordinates_path <- args[[1]]
variance_path <- args[[2]]
output_prefix <- args[[3]]
organelle <- args[[4]]

coordinates <- read.delim(
  coordinates_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)
variance <- read.delim(
  variance_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)
require_columns(
  coordinates,
  c("sample_id", "pc1", "pc2", "species"),
  "PCA coordinate table"
)
require_columns(
  variance,
  c("component", "explained_variance_ratio"),
  "PCA variance table"
)
if (!nrow(coordinates)) {
  stop("PCA coordinate table is empty", call. = FALSE)
}

component_variance <- setNames(
  as.numeric(variance$explained_variance_ratio),
  variance$component
)
if (!all(c("PC1", "PC2") %in% names(component_variance))) {
  stop("PCA variance table must contain PC1 and PC2", call. = FALSE)
}

coordinates$species_group <- normalize_species(coordinates$species)

plot_value <- ggplot(
  coordinates,
  aes(x = as.numeric(pc1), y = as.numeric(pc2), color = species_group)
) +
  geom_hline(yintercept = 0, color = "#888888", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#888888", linewidth = 0.35) +
  geom_point(size = 2.4, alpha = 0.82) +
  scale_color_manual(
    values = species_palette,
    breaks = names(species_palette),
    drop = FALSE,
    name = "Species group"
  ) +
  labs(
    title = paste(organelle, "PCA"),
    subtitle = "Haploid organelle SNP clustering; each point is one sample",
    x = sprintf("PC1 (%.2f%% variance)", component_variance[["PC1"]] * 100),
    y = sprintf("PC2 (%.2f%% variance)", component_variance[["PC2"]] * 100),
    caption = "Colors show species metadata; unresolved metadata is gray."
  ) +
  dudleya_theme(12) +
  theme(legend.position = "right")

save_figure_formats(plot_value, output_prefix, width = 10, height = 7)
