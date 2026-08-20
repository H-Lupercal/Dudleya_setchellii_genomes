suppressPackageStartupMessages(library(ggplot2))

species_palette <- c(
  "D. abramsii ssp. abramsii" = "#0072B2",
  "D. abramsii ssp. bettinae" = "#E69F00",
  "D. abramsii ssp. murina" = "#009E73",
  "D. cymosa" = "#CC79A7",
  "D. setchellii" = "#D55E00",
  "unresolved" = "#777777"
)

normalize_species <- function(values) {
  values[is.na(values) | trimws(values) == ""] <- "unresolved"
  unknown <- setdiff(unique(values), names(species_palette))
  if (length(unknown)) {
    stop(
      paste("no fixed color for species group:", paste(unknown, collapse = ", ")),
      call. = FALSE
    )
  }
  factor(values, levels = names(species_palette))
}

require_columns <- function(value, required, label) {
  missing <- setdiff(required, names(value))
  if (length(missing)) {
    stop(
      paste(label, "is missing columns:", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

dudleya_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "#333333"),
      plot.caption = element_text(color = "#444444", hjust = 0),
      panel.grid.minor = element_blank(),
      legend.key = element_blank()
    )
}

save_figure_formats <- function(plot_value, prefix, width, height) {
  directory <- dirname(prefix)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  png(
    paste0(prefix, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 240,
    bg = "white"
  )
  print(plot_value)
  dev.off()
  pdf(
    paste0(prefix, ".pdf"),
    width = width,
    height = height,
    bg = "white"
  )
  print(plot_value)
  dev.off()
  svg(
    paste0(prefix, ".svg"),
    width = width,
    height = height,
    bg = "white"
  )
  print(plot_value)
  dev.off()
}
