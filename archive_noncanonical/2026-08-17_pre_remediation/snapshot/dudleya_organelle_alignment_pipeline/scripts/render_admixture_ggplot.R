#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "usage: render_admixture_ggplot.R MODE INPUT_TSV OUTPUT_PREFIX ORGANELLE",
    call. = FALSE
  )
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
source(file.path(dirname(normalizePath(script_arg)), "dudleya_plotting.R"))

mode <- args[[1]]
input_path <- args[[2]]
output_prefix <- args[[3]]
organelle <- args[[4]]

cluster_palette <- c(
  "#0072B2", "#E69F00", "#009E73", "#CC79A7",
  "#D55E00", "#56B4E9", "#F0E442", "#6A3D9A"
)

render_structure <- function() {
  q_table <- read.delim(
    input_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )
  require_columns(
    q_table,
    c("sample_id", "plot_group", "popcode"),
    "ADMIXTURE Q table"
  )
  cluster_columns <- grep("^cluster_[0-9]+$", names(q_table), value = TRUE)
  if (!length(cluster_columns)) {
    stop("ADMIXTURE Q table has no cluster columns", call. = FALSE)
  }
  if (!nrow(q_table)) {
    stop("ADMIXTURE Q table is empty", call. = FALSE)
  }

  q_table$sample_index <- seq_len(nrow(q_table))
  long_table <- do.call(
    rbind,
    lapply(seq_along(cluster_columns), function(index) {
      data.frame(
        sample_index = q_table$sample_index,
        cluster = paste0("K", index),
        assignment = as.numeric(q_table[[cluster_columns[[index]]]]),
        stringsAsFactors = FALSE
      )
    })
  )
  long_table$cluster <- factor(
    long_table$cluster,
    levels = paste0("K", seq_along(cluster_columns))
  )

  group_runs <- rle(q_table$plot_group)
  group_ends <- cumsum(group_runs$lengths)
  group_starts <- c(1, head(group_ends, -1) + 1)
  group_midpoints <- (group_starts + group_ends) / 2
  group_labels <- vapply(
    seq_along(group_starts),
    function(index) {
      value <- q_table$popcode[[group_starts[[index]]]]
      if (is.na(value) || value == "") group_runs$values[[index]] else value
    },
    character(1)
  )

  plot_value <- ggplot(long_table, aes(sample_index, assignment, fill = cluster)) +
    geom_col(width = 1, color = NA) +
    geom_vline(
      xintercept = head(group_ends, -1) + 0.5,
      color = "white",
      linewidth = 0.35
    ) +
    scale_fill_manual(
      values = setNames(
        rep(cluster_palette, length.out = length(cluster_columns)),
        paste0("K", seq_along(cluster_columns))
      ),
      name = "Inferred cluster\n(labels arbitrary)"
    ) +
    scale_x_continuous(
      breaks = group_midpoints,
      labels = group_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0)),
      labels = function(values) sprintf("%.1f", values)
    ) +
    labs(
      title = sprintf(
        "%s ADMIXTURE-style organelle clustering, K=%d",
        organelle,
        length(cluster_columns)
      ),
      subtitle = "Samples are ordered by metadata group; white lines separate groups",
      x = "Population code",
      y = "Assignment proportion",
      caption = paste(
        "Cluster numbers and colors are arbitrary labels, not named populations.",
        "Organelle sites form one linked haploid lineage."
      )
    ) +
    dudleya_theme(11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      legend.position = "right"
    )

  save_figure_formats(plot_value, output_prefix, width = 14, height = 7.5)
}

render_cv <- function() {
  summary <- read.delim(
    input_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )
  require_columns(summary, c("organelle", "k", "cv_error"), "ADMIXTURE summary")
  summary <- summary[summary$organelle == organelle, , drop = FALSE]
  if (!nrow(summary)) {
    stop(paste("no ADMIXTURE rows for", organelle), call. = FALSE)
  }
  summary$k <- as.integer(summary$k)
  summary$cv_error <- as.numeric(summary$cv_error)
  split_values <- split(summary$cv_error, summary$k)
  cv_table <- data.frame(
    k = as.integer(names(split_values)),
    mean_cv = vapply(split_values, mean, numeric(1)),
    sd_cv = vapply(
      split_values,
      function(values) if (length(values) > 1) stats::sd(values) else 0,
      numeric(1)
    ),
    replicate_count = vapply(split_values, length, integer(1))
  )
  cv_table <- cv_table[order(cv_table$k), ]
  best_index <- which.min(cv_table$mean_cv)
  best_k <- cv_table$k[[best_index]]
  best_mean <- cv_table$mean_cv[[best_index]]

  plot_value <- ggplot(cv_table, aes(k, mean_cv)) +
    geom_errorbar(
      aes(ymin = pmax(0, mean_cv - sd_cv), ymax = mean_cv + sd_cv),
      width = 0.12,
      color = "#555555"
    ) +
    geom_line(color = "#0072B2", linewidth = 0.8) +
    geom_point(color = "#0072B2", size = 2.7) +
    geom_point(
      data = cv_table[best_index, , drop = FALSE],
      color = "#D55E00",
      size = 4
    ) +
    annotate(
      "label",
      x = best_k,
      y = best_mean,
      label = paste0("Selected K = ", best_k),
      hjust = if (best_k == max(cv_table$k)) 1.1 else -0.1,
      vjust = -0.8,
      size = 3.5
    ) +
    scale_x_continuous(breaks = cv_table$k) +
    labs(
      title = paste(organelle, "ADMIXTURE cross-validation"),
      subtitle = "Lower cross-validation error is better",
      x = "Number of inferred clusters (K)",
      y = "Mean cross-validation error",
      caption = paste0(
        "Points are means across ",
        max(cv_table$replicate_count),
        " replicate(s); error bars show ±1 SD. Orange marks the minimum mean."
      )
    ) +
    dudleya_theme(12)

  save_figure_formats(plot_value, output_prefix, width = 8, height = 6)
}

if (mode == "structure") {
  render_structure()
} else if (mode == "cv") {
  render_cv()
} else {
  stop(paste("unknown ADMIXTURE plot mode:", mode), call. = FALSE)
}
