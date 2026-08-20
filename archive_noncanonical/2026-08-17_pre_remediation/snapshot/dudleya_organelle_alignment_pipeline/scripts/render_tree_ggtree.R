#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    paste(
      "usage: render_tree_ggtree.R TREEFILE METADATA_TSV OUTPUT_PREFIX",
      "ORGANELLE MODE BOOTSTRAP_REPLICATES"
    ),
    call. = FALSE
  )
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
source(file.path(dirname(normalizePath(script_arg)), "dudleya_plotting.R"))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggtree))

tree_path <- args[[1]]
metadata_path <- args[[2]]
output_prefix <- args[[3]]
organelle <- args[[4]]
mode <- args[[5]]
bootstrap_replicates <- as.integer(args[[6]])

if (!mode %in% c("initial", "bootstrap")) {
  stop(paste("unknown tree plot mode:", mode), call. = FALSE)
}

tree <- read.tree(tree_path)
if (is.null(tree) || !length(tree$tip.label)) {
  stop("treefile contains no tips", call. = FALSE)
}
metadata <- read.delim(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)
require_columns(metadata, c("sample_id", "species"), "sample metadata")

tip_species <- metadata$species[match(tree$tip.label, metadata$sample_id)]
tip_metadata <- data.frame(
  label = tree$tip.label,
  species_group = normalize_species(tip_species),
  stringsAsFactors = FALSE
)

tip_count <- length(tree$tip.label)
tip_text_size <- if (tip_count > 150) 1.45 else if (tip_count > 60) 2 else 3
support_text_size <- if (tip_count > 150) 1.25 else 2.3
figure_height <- max(6, min(50, tip_count * 0.13 + 2))

subtitle <- if (mode == "bootstrap") {
  sprintf(
    "Maximum-likelihood tree; internal numbers are UFBoot support percentages (%d replicates)",
    bootstrap_replicates
  )
} else {
  "Maximum-likelihood tree; branch lengths show substitutions per site"
}

plot_value <- ggtree(tree, layout = "rectangular", linewidth = 0.32) %<+%
  tip_metadata
plot_value <- plot_value +
  geom_tiplab(
    aes(color = species_group),
    size = tip_text_size,
    align = FALSE,
    linesize = 0.2,
    show.legend = FALSE
  ) +
  geom_tippoint(
    aes(color = species_group),
    size = 0.01,
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = species_palette,
    breaks = names(species_palette),
    drop = FALSE,
    na.value = species_palette[["unresolved"]],
    name = "Species group"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 15, size = 3, alpha = 1)
    )
  ) +
  geom_treescale(
    x = 0,
    y = -1,
    fontsize = if (tip_count > 150) 2.2 else 3.2,
    linesize = 0.5
  ) +
  labs(
    title = paste(organelle, "phylogenetic tree (ggtree)"),
    subtitle = subtitle,
    x = "Substitutions per site",
    caption = paste(
      "Tip-label colors are metadata annotations and do not affect tree inference.",
      "Branch lengths are substitutions per site."
    )
  ) +
  theme_tree2(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "#333333"),
    plot.caption = element_text(size = 8, color = "#444444", hjust = 0),
    legend.position = "right",
    legend.key = element_blank(),
    axis.title.x = element_text(size = 10)
  )

if (mode == "bootstrap") {
  plot_value <- plot_value +
    geom_text2(
      aes(subset = !isTip & !is.na(label) & label != "", label = label),
      size = support_text_size,
      color = "#333333",
      hjust = 1.15,
      vjust = -0.35
    )
}

plot_value <- plot_value + hexpand(0.35, direction = 1)
save_figure_formats(plot_value, output_prefix, width = 14, height = figure_height)
