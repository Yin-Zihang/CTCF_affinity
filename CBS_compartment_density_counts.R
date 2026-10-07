#!/usr/bin/env Rscript

# Create high-resolution smoothScatter figures from the processed peak tables.
# No ggplot2 or ggpubr installation is required.

read_options <- function(input) {
  if (length(input) == 1L && input[1] == "--version") {
    cat("13review_compartment_density_counts.R 2.0\n")
    quit(save = "no", status = 0)
  }
  if (length(input) %% 2L != 0L || length(input) == 0L) {
    stop("Use --table FILE --output FILE [--layout pages|grid|combined] [--max-chip N] [--max-emsa N] [--nbin N] [--dpi N]")
  }
  keys <- input[seq(1L, length(input), by = 2L)]
  values <- input[seq(2L, length(input), by = 2L)]
  allowed <- c("--table", "--output", "--layout", "--max-chip", "--max-emsa", "--nbin", "--dpi")
  if (any(!keys %in% allowed) || anyDuplicated(keys) || any(!startsWith(keys, "--"))) {
    stop("Unknown or duplicated option: ", paste(keys[!keys %in% allowed], collapse = ", "))
  }
  options <- as.list(setNames(values, substring(keys, 3L)))
  if (is.null(options$table) || is.null(options$output)) {
    stop("Both --table and --output are required")
  }
  options$layout <- if (is.null(options$layout)) "pages" else options$layout
  options$max.emsa <- if (is.null(options$`max-emsa`)) 7 else as.numeric(options$`max-emsa`)
  options$max.chip <- if (is.null(options$`max-chip`)) NA_real_ else as.numeric(options$`max-chip`)
  options$nbin <- if (is.null(options$nbin)) 1600L else as.integer(options$nbin)
  options$dpi <- if (is.null(options$dpi)) 600L else as.integer(options$dpi)
  if (!options$layout %in% c("pages", "grid", "combined")) stop("Invalid --layout")
  if (!is.finite(options$max.emsa) || options$max.emsa <= 0) stop("Invalid --max-emsa")
  if (!is.na(options$max.chip) && (!is.finite(options$max.chip) || options$max.chip <= 0)) stop("Invalid --max-chip")
  if (is.na(options$nbin) || options$nbin < 128L) stop("--nbin must be at least 128")
  if (is.na(options$dpi) || options$dpi < 72L) stop("--dpi must be at least 72")
  options
}

count_ranges <- function(data, label, max_chip, max_emsa) {
  finite <- is.finite(data$ChIP) & is.finite(data$EMSA)
  inside_chip <- finite & data$ChIP >= 0 & data$ChIP <= max_chip
  inside_emsa <- finite & data$EMSA >= 0 & data$EMSA <= max_emsa
  data.frame(category = label, total_peak_rows = nrow(data),
             nonfinite = sum(!finite),
             shown_within_axes = sum(inside_chip & inside_emsa),
             outside_chip_only = sum(finite & !inside_chip & inside_emsa),
             outside_emsa_only = sum(finite & inside_chip & !inside_emsa),
             outside_both = sum(finite & !inside_chip & !inside_emsa))
}

draw_density <- function(data, label, max_chip, max_emsa, nbin) {
  counts <- count_ranges(data, label, max_chip, max_emsa)
  included <- is.finite(data$ChIP) & is.finite(data$EMSA) &
    data$ChIP >= 0 & data$ChIP <= max_chip &
    data$EMSA >= 0 & data$EMSA <= max_emsa
  data <- data[included, , drop = FALSE]
  heading <- sprintf("%s: total %d, shown %d", label,
                     counts$total_peak_rows, counts$shown_within_axes)
  if (nrow(data) < 2L) {
    graphics::plot(NA_real_, NA_real_, type = "n", xlim = c(0, max_chip),
                   ylim = c(0, max_emsa), xlab = "CTCF ChIP RPKM",
                   ylab = "EMSA", main = heading, cex.main = 0.8)
    if (nrow(data) == 1L) graphics::points(data$ChIP, data$EMSA, pch = 20)
    return(invisible(NULL))
  }
  graphics::smoothScatter(data$ChIP, data$EMSA,
                          xlim = c(0, max_chip), ylim = c(0, max_emsa),
                          xlab = "CTCF ChIP RPKM", ylab = "EMSA",
                          main = heading, cex.main = 0.8, nbin = c(nbin, nbin),
                          nrpoints = 100, col = "black", pch = ".",
                          useRaster = TRUE)
  invisible(NULL)
}

main <- function() {
  options <- read_options(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("KernSmooth", quietly = TRUE)) {
    stop("R package 'KernSmooth' is required by smoothScatter")
  }
  table <- utils::read.delim(options$table, check.names = FALSE,
                             stringsAsFactors = FALSE)
  if (!all(c("ChIP", "EMSA") %in% names(table))) {
    stop("The input TSV needs ChIP and EMSA columns")
  }
  table$ChIP <- as.numeric(table$ChIP)
  table$EMSA <- as.numeric(table$EMSA)
  if ("compartments" %in% names(table)) {
    category <- "compartments"
    labels <- c("1.A_Compartment", "2.B_Compartment")
  } else if ("subcompartment" %in% names(table)) {
    category <- "subcompartment"
    labels <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")
  } else {
    stop("The input TSV needs compartments or subcompartment")
  }
  if (nrow(table) == 0L || !any(is.finite(table$ChIP))) {
    stop("No finite ChIP observations found")
  }
  if (is.na(options$max.chip)) {
    options$max.chip <- max(1, table$ChIP[is.finite(table$ChIP)], na.rm = TRUE)
  }
  audit_labels <- c(labels, "All peaks")
  audit <- do.call(rbind, lapply(audit_labels, function(label) {
    selected <- if (label == "All peaks") table else
      table[!is.na(table[[category]]) & table[[category]] == label, , drop = FALSE]
    count_ranges(selected, label, options$max.chip, options$max.emsa)
  }))
  if (options$layout == "combined") labels <- "All peaks"

  output <- options$output
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  extension <- tolower(tools::file_ext(output))
  if (!extension %in% c("pdf", "png")) stop("Output must end in .pdf or .png")
  if (extension == "png" && options$layout == "pages") {
    stop("Use --layout grid or combined for one PNG file; use PDF for multiple pages")
  }
  audit_file <- paste0(tools::file_path_sans_ext(output), "_count_audit.tsv")
  utils::write.table(audit, audit_file, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  print(audit, row.names = FALSE)
  if (options$layout == "grid" && length(labels) == 8L) {
    width <- 16; height <- 9
  } else if (options$layout == "grid" && length(labels) == 2L) {
    width <- 12; height <- 6
  } else {
    width <- 8; height <- 7
  }
  if (extension == "pdf") {
    grDevices::pdf(output, width = width, height = height,
                   onefile = TRUE, useDingbats = FALSE)
  } else if (capabilities("cairo")) {
    grDevices::png(output, width = width, height = height, units = "in",
                   res = options$dpi, type = "cairo", bg = "white")
  } else {
    grDevices::png(output, width = width, height = height, units = "in",
                   res = options$dpi, bg = "white")
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  if (options$layout == "grid" && length(labels) == 8L) {
    graphics::par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
  } else if (options$layout == "grid" && length(labels) == 2L) {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  } else {
    graphics::par(mfrow = c(1, 1), mar = c(5, 5, 3, 2))
  }
  for (label in labels) {
    selected <- if (options$layout == "combined") table else
      table[!is.na(table[[category]]) & table[[category]] == label, , drop = FALSE]
    draw_density(selected, label, options$max.chip, options$max.emsa, options$nbin)
  }
  message("R smoothScatter HD: ", normalizePath(options$table),
          " -> ", file.path(normalizePath(dirname(output)), basename(output)),
          " (nbin=", options$nbin,
          ", layout=", options$layout, ")")
  message("Count audit: ", normalizePath(audit_file))
}

main()
