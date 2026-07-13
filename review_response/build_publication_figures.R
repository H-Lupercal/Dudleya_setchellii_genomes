#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtree)
  library(patchwork)
  library(readr)
  library(stringr)
})

options(warn = 1)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1) {
  stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", script_arg))
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
review_dir <- file.path(repo_root, "review_response")
full_run <- file.path(repo_root, "full_pipeline_run", "results")
output_dir <- file.path(review_dir, "publication_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_order <- c("ABAB", "ABBE", "ABMU", "DUSE", "DUCY", "Other / legacy IDs")
group_colors <- c(
  "ABAB" = "#7B3294",
  "ABBE" = "#E6AB02",
  "ABMU" = "#009E73",
  "DUSE" = "#0072B2",
  "DUCY" = "#D55E00",
  "Other / legacy IDs" = "#777777",
  "Mixed outgroup" = "#3F3F3F"
)
group_shapes <- c(
  "ABAB" = 15,
  "ABBE" = 18,
  "ABMU" = 17,
  "DUSE" = 16,
  "DUCY" = 17,
  "Other / legacy IDs" = 4
)

metadata_path <- file.path(full_run, "07_downstream_sample_set", "included_samples.tsv")
metadata <- read_tsv(metadata_path, show_col_types = FALSE, progress = FALSE) |>
  mutate(
    display_group = case_when(
      str_starts(sample_id, "ABAB_") | str_starts(popcode, "ABAB") ~ "ABAB",
      str_starts(sample_id, "ABBE_") | str_starts(popcode, "ABBE") ~ "ABBE",
      str_starts(sample_id, "ABMU_") | str_starts(popcode, "ABMU") ~ "ABMU",
      str_detect(species, fixed("setchellii")) ~ "DUSE",
      str_detect(species, fixed("cymosa")) ~ "DUCY",
      TRUE ~ "Other / legacy IDs"
    ),
    display_group = factor(display_group, levels = group_order),
    popcode = na_if(popcode, "")
  )

if (nrow(metadata) != 278 || anyDuplicated(metadata$sample_id)) {
  stop("Expected exactly 278 unique included samples.")
}

tree_paths <- c(
  cpDNA = file.path(review_dir, "cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.treefile"),
  mtDNA = file.path(review_dir, "mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.treefile")
)
trees <- lapply(tree_paths, read.tree)

for (organelle in names(trees)) {
  tree <- trees[[organelle]]
  if (Ntip(tree) != 278) {
    stop(organelle, " tree has ", Ntip(tree), " tips; expected 278.")
  }
  if (!setequal(tree$tip.label, metadata$sample_id)) {
    stop(organelle, " tree tips do not match included_samples.tsv.")
  }
  if (length(tree$node.label) != tree$Nnode) {
    stop(organelle, " tree does not contain one support field per internal node.")
  }
}

if (!setequal(trees$cpDNA$tip.label, trees$mtDNA$tip.label)) {
  stop("cpDNA and mtDNA trees do not contain the same tips.")
}

descendant_tip_indices <- function(tree, node) {
  tips <- integer()
  queue <- node
  while (length(queue) > 0) {
    current <- queue[[1]]
    queue <- queue[-1]
    children <- tree$edge[tree$edge[, 1] == current, 2]
    tips <- c(tips, children[children <= Ntip(tree)])
    queue <- c(queue, children[children > Ntip(tree)])
  }
  sort(unique(tips))
}

monophyletic_node <- function(tree, tip_names) {
  tip_names <- intersect(tip_names, tree$tip.label)
  if (length(tip_names) < 2) {
    return(NA_integer_)
  }
  node <- getMRCA(tree, tip_names)
  if (is.null(node) || is.na(node)) {
    return(NA_integer_)
  }
  descendants <- tree$tip.label[descendant_tip_indices(tree, node)]
  if (!setequal(descendants, tip_names)) {
    return(NA_integer_)
  }
  as.integer(node)
}

collapse_candidates <- function(tree, metadata) {
  candidates <- list()
  outgroup_tips <- metadata |>
    filter(as.character(display_group) %in% c("ABAB", "ABMU")) |>
    pull(sample_id)
  outgroup_node <- monophyletic_node(tree, outgroup_tips)

  if (!is.na(outgroup_node)) {
    candidates[[length(candidates) + 1]] <- tibble(
      node = outgroup_node,
      clade = "ABAB + ABMU outgroup",
      sample_count = length(outgroup_tips),
      display_group = "Mixed outgroup"
    )
  }

  population_groups <- metadata |>
    filter(!is.na(popcode), !(sample_id %in% outgroup_tips)) |>
    group_by(popcode, display_group) |>
    summarise(sample_ids = list(sample_id), sample_count = n(), .groups = "drop") |>
    filter(sample_count >= 3)

  for (i in seq_len(nrow(population_groups))) {
    row <- population_groups[i, ]
    node <- monophyletic_node(tree, row$sample_ids[[1]])
    if (!is.na(node)) {
      candidates[[length(candidates) + 1]] <- tibble(
        node = node,
        clade = str_replace_all(row$popcode, "_", " "),
        sample_count = row$sample_count,
        display_group = as.character(row$display_group)
      )
    }
  }

  if (length(candidates) == 0) {
    return(tibble(node = integer(), clade = character(), sample_count = integer(), display_group = character()))
  }
  bind_rows(candidates) |>
    distinct(node, .keep_all = TRUE) |>
    arrange(desc(sample_count), clade)
}

node_annotations <- function(tree) {
  tibble(
    node = seq.int(Ntip(tree) + 1, Ntip(tree) + tree$Nnode),
    bootstrap = suppressWarnings(as.numeric(tree$node.label))
  )
}

make_collapsed_tree <- function(tree, organelle, metadata) {
  p <- ggtree(tree, size = 0.32, color = "#333333")
  p$data <- p$data |>
    left_join(
      metadata |> select(sample_id, display_group, popcode),
      by = c("label" = "sample_id")
    ) |>
    left_join(node_annotations(tree), by = "node") |>
    mutate(
      display_group = as.character(display_group),
      collapsed_label = NA_character_,
      collapsed_group = NA_character_
    )

  collapsed <- collapse_candidates(tree, metadata)
  for (i in seq_len(nrow(collapsed))) {
    row <- collapsed[i, ]
    fill_color <- unname(group_colors[[row$display_group]])
    p <- collapse(
      p,
      node = row$node,
      mode = "max",
      fill = fill_color,
      color = fill_color,
      alpha = 0.68,
      size = 0.25
    )
    p$data$collapsed_label[p$data$node == row$node] <- paste0(
      row$clade, " (n=", row$sample_count, ")"
    )
    p$data$collapsed_group[p$data$node == row$node] <- row$display_group
  }

  max_x <- max(p$data$x, na.rm = TRUE)
  outgroup_points <- p$data |>
    filter(isTip, display_group %in% c("ABAB", "ABMU"), !is.na(x), !is.na(y))
  outgroup_annotation <- tibble(
    x = max(outgroup_points$x) + max_x * 0.018,
    y = mean(range(outgroup_points$y)),
    label = paste0("ABAB + ABMU outgroup\n(n=", nrow(outgroup_points), ")")
  )
  p <- p +
    geom_tippoint(
      aes(color = display_group),
      size = 1.05,
      alpha = 0.88,
      show.legend = TRUE,
      na.rm = TRUE
    ) +
    geom_label(
      data = outgroup_annotation,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.0,
      lineheight = 0.95,
      color = "#333333",
      fill = "white",
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    geom_text2(
      aes(
        subset = !isTip & !is.na(collapsed_label),
        label = collapsed_label,
        color = collapsed_group
      ),
      hjust = -0.04,
      size = 2.25,
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_point2(
      aes(
        subset = !isTip & is.na(collapsed_label) & !is.na(bootstrap) & bootstrap >= 95
      ),
      shape = 21,
      size = 1.15,
      color = "#111111",
      fill = "#111111",
      stroke = 0.15,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    scale_color_manual(values = group_colors, breaks = group_order, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.34))) +
    labs(
      title = organelle,
      subtitle = paste0(
        "Maximum likelihood; rooted with ABAB + ABMU | ",
        "outgroup annotated; black nodes have UFBoot >=95"
      ),
      x = "Substitutions per site",
      color = "Sample group"
    ) +
    theme_tree2() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 7.2, color = "#444444"),
      axis.title.x = element_text(size = 7.5),
      axis.text.x = element_text(size = 6.5),
      legend.position = "bottom",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      plot.margin = margin(5, 10, 5, 5)
    )

  list(plot = p, collapsed = collapsed)
}

make_full_tree <- function(tree, organelle, metadata) {
  p <- ggtree(tree, size = 0.24, color = "#333333")
  p$data <- p$data |>
    left_join(metadata |> select(sample_id, display_group), by = c("label" = "sample_id")) |>
    left_join(node_annotations(tree), by = "node") |>
    mutate(display_group = as.character(display_group))
  max_x <- max(p$data$x, na.rm = TRUE)

  p +
    geom_tiplab(
      aes(color = display_group),
      size = 1.55,
      offset = max_x * 0.008,
      linesize = 0.16
    ) +
    geom_text2(
      aes(
        subset = !isTip & !is.na(bootstrap) & bootstrap >= 95,
        label = sprintf("%.0f", bootstrap)
      ),
      size = 1.25,
      color = "#111111",
      nudge_y = 0.18,
      check_overlap = TRUE
    ) +
    scale_color_manual(values = group_colors, breaks = group_order, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.28))) +
    labs(
      title = paste0(organelle, " rooted maximum-likelihood tree"),
      subtitle = "All 278 samples; ABAB + ABMU outgroup; UFBoot values >=95 shown",
      x = "Substitutions per site",
      color = "Sample group"
    ) +
    theme_tree2() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 7.5, color = "#444444"),
      axis.title.x = element_text(size = 7.5),
      axis.text.x = element_text(size = 6.5),
      legend.position = "bottom",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.5),
      plot.margin = margin(5, 12, 5, 5)
    )
}

read_variance <- function(organelle) {
  path <- file.path(full_run, "15_pca", paste0(organelle, ".primary.pca.variance.tsv"))
  values <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  setNames(values$explained_variance_ratio, values$component)
}

make_pca_panel <- function(organelle, snp_count) {
  path <- file.path(review_dir, paste0(organelle, ".primary.pca.requested_groups.coordinates.tsv"))
  coords <- read_tsv(path, show_col_types = FALSE, progress = FALSE) |>
    mutate(
      display_group = factor(display_group, levels = group_order),
      popcode = na_if(popcode, "")
    )
  if (nrow(coords) != 278 || anyDuplicated(coords$sample_id)) {
    stop(organelle, " PCA coordinates do not contain 278 unique samples.")
  }

  centroids <- coords |>
    filter(!is.na(popcode)) |>
    group_by(popcode, display_group) |>
    summarise(pc1 = mean(pc1), pc2 = mean(pc2), sample_count = n(), .groups = "drop") |>
    filter(sample_count >= 2) |>
    mutate(label = str_replace_all(popcode, "_", " "))

  selected_labels <- coords |>
    filter(sample_id %in% c("ABAB_MAD_LP_222_Du-589", "DU-443"))
  variance <- read_variance(organelle)

  ggplot(coords, aes(pc1, pc2)) +
    geom_hline(yintercept = 0, color = "#B0B0B0", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "#B0B0B0", linewidth = 0.3) +
    geom_point(
      aes(color = display_group, shape = display_group),
      size = 1.75,
      alpha = 0.72,
      stroke = 0.25
    ) +
    geom_point(
      data = centroids,
      aes(fill = display_group),
      shape = 21,
      size = 2.5,
      color = "white",
      stroke = 0.5,
      show.legend = FALSE
    ) +
    geom_text_repel(
      data = centroids,
      aes(label = label, color = display_group),
      size = 2.0,
      box.padding = 0.22,
      point.padding = 0.12,
      min.segment.length = 0.05,
      segment.size = 0.2,
      max.overlaps = Inf,
      seed = 42,
      force = 0.6,
      show.legend = FALSE
    ) +
    geom_label_repel(
      data = selected_labels,
      aes(label = sample_id, color = display_group),
      size = 1.9,
      fill = "white",
      label.size = 0.15,
      min.segment.length = 0,
      seed = 42,
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_colors, breaks = group_order, drop = FALSE) +
    scale_fill_manual(values = group_colors, breaks = group_order, drop = FALSE) +
    scale_shape_manual(values = group_shapes, breaks = group_order, drop = FALSE) +
    labs(
      title = organelle,
      subtitle = paste0("278 samples | ", format(snp_count, big.mark = ","), " SNPs"),
      x = sprintf("PC1 (%.2f%%)", variance[["PC1"]] * 100),
      y = sprintf("PC2 (%.2f%%)", variance[["PC2"]] * 100),
      color = "Group",
      shape = "Group"
    ) +
    theme_classic(base_size = 9) +
    theme(
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.45),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 7.5, color = "#555555"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      legend.position = "bottom",
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7),
      plot.margin = margin(5, 8, 5, 5)
    )
}

make_branch_comparison <- function() {
  shared_path <- file.path(review_dir, "strongly_supported_shared_splits.tsv")
  conflicts_path <- file.path(review_dir, "strongly_supported_conflicting_splits.tsv")
  shared <- read_tsv(shared_path, show_col_types = FALSE, progress = FALSE)
  conflicts <- read_tsv(conflicts_path, show_col_types = FALSE, progress = FALSE)

  shared_labels <- c(
    "ABBE + CY BAL/BGL (n=9)",
    "CY ALA (n=3)",
    "CY BAL (n=3)",
    "COM (n=9)",
    "SILG (n=4)",
    "CY HICN/SIE (n=9)",
    "DOV (n=17)",
    "ABAB + ABMU (n=13)"
  )
  shared_display <- shared |>
    slice_head(n = length(shared_labels)) |>
    mutate(branch = factor(shared_labels, levels = rev(shared_labels)))

  shared_plot <- ggplot(shared_display, aes(y = branch)) +
    geom_segment(
      aes(x = cpDNA_support, xend = mtDNA_support, yend = branch),
      color = "#9E9E9E",
      linewidth = 0.6
    ) +
    geom_point(aes(x = cpDNA_support, color = "cpDNA"), size = 3.0) +
    geom_point(aes(x = mtDNA_support, color = "mtDNA"), size = 3.0, shape = 17) +
    geom_text(aes(x = cpDNA_support, label = cpDNA_support), color = "white", size = 2.1) +
    geom_text(aes(x = mtDNA_support, label = mtDNA_support), color = "white", size = 2.1) +
    scale_color_manual(values = c("cpDNA" = "#0072B2", "mtDNA" = "#D55E00")) +
    scale_x_continuous(limits = c(94, 101), breaks = 95:100, expand = expansion(mult = c(0.02, 0.03))) +
    labs(
      title = "A  Strongly supported branches shared by both organelles",
      x = "Ultrafast bootstrap support (%)",
      y = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 7.5),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      legend.position = "top",
      legend.text = element_text(size = 7.5),
      plot.margin = margin(5, 8, 5, 5)
    )

  conflict_labels <- tibble(
    row_id = seq_len(min(3, nrow(conflicts))),
    cp_label = c(
      "TIL pair (n=2)",
      "DAN-centered clade (n=12)",
      "HICN/SIE composition A (n=8)"
    )[seq_len(min(3, nrow(conflicts)))],
    mt_label = c(
      "QUI2-centered DUSE clade (n=22)",
      "SILR clade (n=6)",
      "HICN/SIE composition B (n=8)"
    )[seq_len(min(3, nrow(conflicts)))]
  ) |>
    bind_cols(conflicts |> slice_head(n = 3)) |>
    mutate(
      y = rev(seq_len(n())),
      cp_text = paste0(cp_label, "\nUFBoot ", source_support),
      mt_text = paste0(mt_label, "\nUFBoot ", conflicting_support)
    )

  conflict_plot <- ggplot(conflict_labels) +
    geom_segment(
      aes(x = 0.28, xend = 0.72, y = y, yend = y),
      color = "#B2182B",
      linewidth = 0.7,
      linetype = "22"
    ) +
    geom_label(
      aes(x = 0.20, y = y, label = cp_text),
      fill = "#E6F2F8",
      color = "#004C6D",
      label.size = 0.25,
      size = 2.65,
      lineheight = 0.95
    ) +
    geom_label(
      aes(x = 0.80, y = y, label = mt_text),
      fill = "#FBE9E1",
      color = "#8A2D0A",
      label.size = 0.25,
      size = 2.65,
      lineheight = 0.95
    ) +
    annotate("text", x = 0.20, y = max(conflict_labels$y) + 0.65, label = "cpDNA split", fontface = "bold", size = 3.1) +
    annotate("text", x = 0.80, y = max(conflict_labels$y) + 0.65, label = "Incompatible mtDNA split", fontface = "bold", size = 3.1) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.4, max(conflict_labels$y) + 0.9), clip = "off") +
    labs(
      title = "B  Strongly supported incompatible branches",
      subtitle = "Each pair represents two bipartitions that cannot both occur in one tree",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7.5, color = "#555555"),
      plot.margin = margin(5, 15, 5, 15)
    )

  shared_plot / conflict_plot +
    plot_layout(heights = c(1.35, 1)) +
    plot_annotation(
      title = "cpDNA and mtDNA branch concordance",
      subtitle = "Stage 19 maximum-likelihood trees; strong support defined as UFBoot >=95",
      caption = "The mtDNA alignment has 146 informative sites; the cpDNA alignment has 2,022.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.5, color = "#444444"),
        plot.caption = element_text(size = 7, color = "#555555")
      )
    )
}

save_plot <- function(plot, stem, width, height, dpi = 600) {
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(output_dir, paste0(stem, ".svg"))
  png_path <- file.path(output_dir, paste0(stem, ".png"))

  cairo_pdf(pdf_path, width = width, height = height, family = "sans")
  print(plot)
  dev.off()

  svg(svg_path, width = width, height = height, family = "sans", onefile = FALSE)
  print(plot)
  dev.off()

  png(
    png_path,
    width = round(width * dpi),
    height = round(height * dpi),
    res = dpi,
    type = "cairo",
    antialias = "subpixel"
  )
  print(plot)
  dev.off()

  tibble(
    figure = stem,
    pdf = basename(pdf_path),
    svg = basename(svg_path),
    png = basename(png_path),
    width_in = width,
    height_in = height,
    dpi = dpi
  )
}

collapsed_cp <- make_collapsed_tree(trees$cpDNA, "cpDNA", metadata)
collapsed_mt <- make_collapsed_tree(trees$mtDNA, "mtDNA", metadata)

tree_plate <- collapsed_cp$plot + collapsed_mt$plot +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Rooted organelle phylogenies",
    subtitle = "ABAB + ABMU outgroup; GTR+F+G4; 1,000 ultrafast bootstrap replicates with BNNI",
    caption = "Population clades are collapsed only when monophyletic. Full tip labels and numerical supports are in the supplementary trees.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5, color = "#444444"),
      plot.caption = element_text(size = 7, color = "#555555")
    )
  ) & theme(legend.position = "bottom")

pca_plate <- make_pca_panel("cpDNA", 2022) + make_pca_panel("mtDNA", 146) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Organelle genetic structure",
    subtitle = "Points are samples; labeled white-centered points are population centroids",
    caption = "PCA is based on filtered haploid SNP alignments. Missing states were mean-imputed per retained site.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5, color = "#444444"),
      plot.caption = element_text(size = 7, color = "#555555")
    )
  ) & theme(legend.position = "bottom")

branch_comparison <- make_branch_comparison()
full_cp <- make_full_tree(trees$cpDNA, "cpDNA", metadata)
full_mt <- make_full_tree(trees$mtDNA, "mtDNA", metadata)

manifest <- bind_rows(
  save_plot(tree_plate, "figure_1_rooted_collapsed_trees", 14, 9),
  save_plot(branch_comparison, "figure_2_cpdna_mtdna_branch_comparison", 9, 8),
  save_plot(pca_plate, "figure_3_pca_requested_groups", 13, 6.8),
  save_plot(full_cp, "supplementary_figure_1_cpDNA_full_tree", 8.5, 24),
  save_plot(full_mt, "supplementary_figure_2_mtDNA_full_tree", 8.5, 24)
)

collapsed_summary <- bind_rows(
  collapsed_cp$collapsed |> mutate(organelle = "cpDNA", .before = 1),
  collapsed_mt$collapsed |> mutate(organelle = "mtDNA", .before = 1)
)

write_tsv(manifest, file.path(output_dir, "figure_manifest.tsv"))
write_tsv(collapsed_summary, file.path(output_dir, "collapsed_clades.tsv"))

expected <- unlist(manifest[, c("pdf", "svg", "png")], use.names = FALSE)
expected_paths <- file.path(output_dir, expected)
if (any(!file.exists(expected_paths)) || any(file.info(expected_paths)$size <= 0)) {
  stop("One or more publication figure files are missing or empty.")
}

message("Publication figures written to: ", output_dir)
message("Collapsed cpDNA clades: ", nrow(collapsed_cp$collapsed))
message("Collapsed mtDNA clades: ", nrow(collapsed_mt$collapsed))
