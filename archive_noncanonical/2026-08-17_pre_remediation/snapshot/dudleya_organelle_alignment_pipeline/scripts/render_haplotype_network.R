#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "usage: render_haplotype_network.R FASTA METADATA PREFIX ORGANELLE",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(pegas))

fasta_path <- args[[1]]
metadata_path <- args[[2]]
prefix <- args[[3]]
organelle <- args[[4]]

species_palette <- c(
  "D. abramsii ssp. abramsii" = "#0072B2",
  "D. abramsii ssp. bettinae" = "#E69F00",
  "D. abramsii ssp. murina" = "#009E73",
  "D. cymosa" = "#CC79A7",
  "D. setchellii" = "#D55E00",
  "unresolved" = "#999999"
)

write_tsv <- function(value, suffix) {
  write.table(
    value,
    file = paste0(prefix, suffix),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

dna <- read.dna(fasta_path, format = "fasta")
metadata <- read.delim(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)
required_metadata <- c("sample_id", "species_group", "popcode")
missing_columns <- setdiff(required_metadata, names(metadata))
if (length(missing_columns)) {
  stop(
    paste("metadata is missing columns:", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}
if (!identical(rownames(dna), metadata$sample_id)) {
  stop("FASTA and metadata sample IDs or order differ", call. = FALSE)
}
metadata$species_group[metadata$species_group == ""] <- "unresolved"
unknown_species <- setdiff(unique(metadata$species_group), names(species_palette))
if (length(unknown_species)) {
  stop(
    paste("no fixed color for species group:", paste(unknown_species, collapse = ", ")),
    call. = FALSE
  )
}

haps <- haplotype(dna)
hap_indexes <- attr(haps, "index")
if (length(hap_indexes) < 2) {
  stop("fewer than two haplotypes supplied to renderer", call. = FALSE)
}
net <- haploNet(haps)
hap_width <- max(3, nchar(as.character(length(hap_indexes))))
haplotype_ids <- sprintf(paste0("H%0", hap_width, "d"), seq_along(hap_indexes))
attr(net, "labels") <- haplotype_ids

sample_haplotype <- rep(NA_character_, nrow(dna))
for (index in seq_along(hap_indexes)) {
  sample_haplotype[hap_indexes[[index]]] <- haplotype_ids[[index]]
}
if (anyNA(sample_haplotype)) {
  stop("could not assign every sample to a haplotype", call. = FALSE)
}

assignments <- data.frame(
  sample_id = metadata$sample_id,
  organelle = organelle,
  haplotype_id = sample_haplotype,
  species_group = metadata$species_group,
  popcode = metadata$popcode,
  stringsAsFactors = FALSE
)

haplotype_summary <- do.call(
  rbind,
  lapply(seq_along(haplotype_ids), function(index) {
    members <- hap_indexes[[index]]
    nonblank_popcodes <- metadata$popcode[members]
    nonblank_popcodes <- nonblank_popcodes[nonblank_popcodes != ""]
    data.frame(
      organelle = organelle,
      haplotype_id = haplotype_ids[[index]],
      sample_count = length(members),
      species_group_count = length(unique(metadata$species_group[members])),
      popcode_count = length(unique(nonblank_popcodes)),
      stringsAsFactors = FALSE
    )
  })
)

edge_rows <- function(edge_matrix, alternative_link) {
  if (is.null(edge_matrix) || !nrow(edge_matrix)) {
    return(NULL)
  }
  data.frame(
    organelle = organelle,
    from_haplotype = haplotype_ids[as.integer(edge_matrix[, 1])],
    to_haplotype = haplotype_ids[as.integer(edge_matrix[, 2])],
    mutation_steps = as.numeric(edge_matrix[, 3]),
    alternative_link = alternative_link,
    stringsAsFactors = FALSE
  )
}
primary_edges <- edge_rows(unclass(net), FALSE)
alternative_edges <- edge_rows(attr(net, "alter.links"), TRUE)
edges <- rbind(primary_edges, alternative_edges)

species_levels <- names(species_palette)
pie_counts <- t(vapply(
  hap_indexes,
  function(members) {
    as.integer(table(factor(
      metadata$species_group[members],
      levels = species_levels
    )))
  },
  integer(length(species_levels))
))
colnames(pie_counts) <- species_levels
rownames(pie_counts) <- haplotype_ids
node_size_scale <- max(1, median(net[, 3]))
node_sizes <- sqrt(vapply(hap_indexes, length, integer(1))) * node_size_scale
observed_species <- species_levels[colSums(pie_counts) > 0]
show_all_labels <- length(haplotype_ids) <= 25
long_edge_threshold <- max(
  5,
  as.numeric(quantile(net[, 3], probs = 0.9, names = FALSE, type = 1))
)

draw_network <- function(coordinates = NULL) {
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  par(mar = c(3, 2, 4, 1), xpd = FALSE)
  plotted <- plot(
    net,
    size = node_sizes,
    pie = pie_counts,
    bg = unname(species_palette),
    labels = show_all_labels,
    legend = FALSE,
    show.mutation = if (show_all_labels) 3 else 0,
    threshold = 0,
    xy = coordinates,
    fast = FALSE
  )
  if (!show_all_labels) {
    long_edges <- which(net[, 3] >= long_edge_threshold)
    if (length(long_edges)) {
      from <- as.integer(net[long_edges, 1])
      to <- as.integer(net[long_edges, 2])
      text(
        (plotted$xx[from] + plotted$xx[to]) / 2,
        (plotted$yy[from] + plotted$yy[to]) / 2,
        labels = net[long_edges, 3],
        cex = 0.65,
        font = 2
      )
    }
  }
  title(main = paste(organelle, "haplotype network (pegas)"))
  mtext(
    if (show_all_labels) {
      "Node area = sample count; sectors = species group; edge labels = mutation steps"
    } else {
      paste0(
        "Node area = sample count; sectors = species group; ",
        "primary-edge labels shown for ≥", long_edge_threshold, " mutation steps"
      )
    },
    side = 1,
    line = 1,
    cex = 0.75
  )
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = observed_species,
    fill = unname(species_palette[observed_species]),
    border = NA,
    bty = "n",
    cex = 0.8,
    title = "Species group"
  )
  invisible(list(x = plotted$xx, y = plotted$yy))
}

set.seed(20260714)
png(
  paste0(prefix, ".haplotype_network.png"),
  width = 2400,
  height = 1800,
  res = 240,
  bg = "white"
)
layout <- draw_network()
dev.off()

pdf(
  paste0(prefix, ".haplotype_network.pdf"),
  width = 10,
  height = 7.5,
  bg = "white"
)
draw_network(layout)
dev.off()

svg(
  paste0(prefix, ".haplotype_network.svg"),
  width = 10,
  height = 7.5,
  bg = "white"
)
draw_network(layout)
dev.off()

layout_table <- data.frame(
  organelle = organelle,
  haplotype_id = haplotype_ids,
  x = layout$x,
  y = layout$y,
  stringsAsFactors = FALSE
)
renderer_summary <- data.frame(
  organelle = organelle,
  sample_count = nrow(dna),
  haplotype_count = length(haplotype_ids),
  edge_count = nrow(edges),
  species_group_count = length(unique(metadata$species_group)),
  stringsAsFactors = FALSE
)

write_tsv(assignments, ".haplotype_assignments.tsv")
write_tsv(haplotype_summary, ".haplotype_summary.tsv")
write_tsv(edges, ".haplotype_network_edges.tsv")
write_tsv(layout_table, ".haplotype_network_layout.tsv")
write_tsv(
  renderer_summary,
  ".haplotype_network_renderer_summary.tsv"
)
