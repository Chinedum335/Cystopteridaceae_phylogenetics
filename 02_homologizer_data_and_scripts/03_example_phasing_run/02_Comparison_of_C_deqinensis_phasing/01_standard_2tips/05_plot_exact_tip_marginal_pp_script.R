# =============================================================================
# Exact-tip marginal posterior-probability plot for Homologizer
# =============================================================================
#
# PURPOSE
# -------
# This script gives every constrained tip its own heatmap cell and reports the
# marginal posterior probability (PP) that the displayed copy occupies that
# exact tip. It never combines A1/A2, B1/B2, or other paired tips.
#
# Use it for:
#   * standard Homologizer analyses with one tip per subgenome; or
#   * paired-tip analyses when exact placement within a pair is of interest.
#
# HOW TO RUN IN RSTUDIO
# ---------------------
# 1. Edit only the USER SETTINGS section below.
# 2. Select the entire script (Ctrl+A).
# 3. Run the selection (Ctrl+Enter), or click Source.
#
# REQUIRED INPUTS
# ---------------
# The analysis directory must contain, or be immediately above, a Homologizer
# output folder containing homologizer_locus_*_phase.log and a map tree. A
# A gene-copy mapping CSV (normally subgenomes_names.csv) with Sample,
# Subgenome, and locus columns must be located in the analysis directory or one
# of its parent directories.

# =============================================================================
# 1. USER SETTINGS -- edit these values before running
# =============================================================================

# Folder containing output1/output2, or the output folder itself. Forward
# slashes work on Windows, macOS, and Linux. "." means the R working directory.
analysis_dir <- "."
# Example:
# analysis_dir <- "C:/Users/your_name/path/to/temp_outputs"

# Use "auto" when there is exactly one Homologizer output folder. Otherwise,
# enter its name (for example, "output2") or its full path.
output_folder <- "auto"

# Destination for the PDF, PNG, and two CSV summaries. "auto" creates an
# exact_tip_plot_results folder inside analysis_dir. You may instead supply a
# full path. A relative path is interpreted relative to analysis_dir.
figure_dir <- "auto"
# Example:
# figure_dir <- "C:/Users/your_name/path/to/figure_results"

# One or more viridis-family palettes. Multiple palettes produce multiple
# figures, for example c("cividis", "plasma", "turbo").
palette_names <- c("plasma")

# Fraction of MCMC samples discarded from the beginning of each phase log.
burnin <- 0.20

# Uniform multiplier applied to every branch length for plotting only.
# 1.00 = original branch lengths; 0.75 = moderately compact; 0.50 = half length.
# The input tree file is never modified.
branch_length_scale <- 1.00

# Empty/unreported assignments are white boxes. These settings control their
# outlines but have no effect on posterior probabilities.
blank_box_border_color <- "grey72"
blank_box_border_width <- 0.30

# =============================================================================
# 2. CHECK SETTINGS AND CREATE THE OUTPUT DIRECTORY
# =============================================================================

analysis_dir <- normalizePath(analysis_dir, winslash = "/", mustWork = TRUE)

# Resolve the destination now that analysis_dir is absolute. This prevents the
# output location from depending on RStudio's working directory during saving.
if (!is.character(figure_dir) || length(figure_dir) != 1 || !nzchar(figure_dir)) {
  stop("figure_dir must be \"auto\" or one non-empty folder path.")
}
if (identical(tolower(figure_dir), "auto")) {
  figure_dir <- file.path(analysis_dir, "exact_tip_plot_results")
} else {
  figure_dir <- path.expand(figure_dir)
  figure_dir_is_absolute <- grepl(
    "^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", figure_dir
  )
  if (!figure_dir_is_absolute) {
    figure_dir <- file.path(analysis_dir, figure_dir)
  }
}

allowed_palettes <- c(
  "viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo"
)
if (any(!palette_names %in% allowed_palettes)) {
  stop("Unknown palette. Choose from: ", paste(allowed_palettes, collapse = ", "))
}
if (!is.finite(burnin) || burnin < 0 || burnin >= 1) {
  stop("burnin must be a number greater than or equal to 0 and less than 1.")
}
if (!is.finite(branch_length_scale) || branch_length_scale <= 0) {
  stop("branch_length_scale must be a number greater than 0.")
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(figure_dir)) {
  stop(
    "R could not create the output directory:\n", figure_dir,
    "\nChoose a location where you have permission to create files."
  )
}
figure_dir <- normalizePath(figure_dir, winslash = "/", mustWork = TRUE)

# Test write access immediately, before the time-consuming plot is assembled.
write_test <- tempfile(pattern = "homologizer_write_test_", tmpdir = figure_dir)
write_test_ok <- tryCatch(file.create(write_test), warning = function(w) FALSE,
                          error = function(e) FALSE)
if (!isTRUE(write_test_ok)) {
  stop(
    "The output directory exists, but R cannot write to it:\n", figure_dir,
    "\nChoose another figure_dir or check the folder permissions."
  )
}
unlink(write_test)
message("Verified writable figure destination: ", figure_dir)

# =============================================================================
# 3. LOCATE HOMOLOGIZER INPUT FILES
# =============================================================================
# The following code finds the phase logs, gene-copy table, and MAP tree. No
# paths below this point normally need to be edited.

# Accept either the analysis directory or the Homologizer output directory.
if (identical(tolower(output_folder), "auto")) {
  candidate_dirs <- unique(c(
    analysis_dir,
    list.dirs(analysis_dir, recursive = FALSE, full.names = TRUE)
  ))
  candidate_dirs <- candidate_dirs[vapply(candidate_dirs, function(path) {
    length(list.files(
      path, pattern = "^homologizer_locus_[0-9]+_phase[.]log$"
    )) > 0
  }, logical(1))]
  if (length(candidate_dirs) != 1) {
    stop(
      "Auto-detection found ", length(candidate_dirs),
      " Homologizer output folders. Supply the folder name explicitly."
    )
  }
  output_dir <- candidate_dirs[1]
} else {
  output_dir <- if (grepl("^([A-Za-z]:)?[/\\\\]", output_folder)) {
    output_folder
  } else {
    file.path(analysis_dir, output_folder)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
}
output_label <- basename(output_dir)
prefix <- file.path(output_dir, "homologizer")

# Find a gene-copy table by its required columns. Searching several ancestors
# lets the script run from temp_outputs, output1/output2, or a copied folder.
search_roots <- unique(c(
  analysis_dir,
  dirname(analysis_dir),
  dirname(dirname(analysis_dir)),
  dirname(dirname(dirname(analysis_dir)))
))
csv_candidates <- unique(c(
  file.path(search_roots, "subgenomes_names.csv"),
  file.path(search_roots, "subgenomes_name.csv"),
  file.path(search_roots, "cystopteridaceae_genomes.csv"),
  unlist(lapply(search_roots, function(path) {
    list.files(path, pattern = "[.]csv$", full.names = TRUE)
  }), use.names = FALSE)
))
csv_candidates <- csv_candidates[file.exists(csv_candidates)]
valid_csv <- csv_candidates[vapply(csv_candidates, function(path) {
  header <- tryCatch(
    names(read.csv(path, nrows = 1, check.names = FALSE)),
    error = function(e) character()
  )
  all(c("Sample", "Subgenome") %in% header)
}, logical(1))]
if (length(valid_csv) == 0) {
  stop(
    "Could not find a gene-copy mapping CSV containing Sample and Subgenome ",
    "columns (normally subgenomes_names.csv)."
  )
}
genecopy_file <- valid_csv[1]

tree_candidates <- c(
  file.path(output_dir, "homologizer_map1.tree"),
  file.path(output_dir, "homologizer_map.tree"),
  list.files(output_dir, pattern = "^homologizer_map.*[.]tree$", full.names = TRUE)
)
tree_candidates <- unique(tree_candidates[file.exists(tree_candidates)])
if (length(tree_candidates) == 0) {
  stop("Could not find homologizer_map1.tree or homologizer_map.tree in ", output_dir)
}
tree_file <- tree_candidates[1]

# =============================================================================
# 4. LOAD PACKAGES AND READ THE ANALYSIS DESIGN
# =============================================================================
# Install missing packages once from the R console before running this script.
# ggtree and treeio are Bioconductor packages; the remaining packages are on
# CRAN. The script stops with an explicit list if anything is unavailable.

required_packages <- c(
  "ggplot2", "magrittr", "tidyr", "dplyr", "ggtree", "treeio", "ape"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop(
    "Install these required R packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(magrittr)
  library(tidyr)
  library(dplyr)
  library(ggtree)
  library(ape)
})

genecopy_map <- read.csv(
  genecopy_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)
if (ncol(genecopy_map) < 3) {
  stop("The gene-copy CSV must have at least one locus column after Sample and Subgenome.")
}
samples <- split(as.character(genecopy_map$Subgenome), genecopy_map$Sample)
loci <- names(genecopy_map)[3:ncol(genecopy_map)]

first_log <- paste0(prefix, "_locus_1_phase.log")
if (!file.exists(first_log)) {
  stop("Missing expected log file: ", first_log)
}
log_columns <- names(read.delim(first_log, nrows = 1, check.names = FALSE))
samples <- samples[vapply(samples, function(tips) {
  length(tips) > 0 && all(tips %in% log_columns)
}, logical(1))]
if (length(samples) == 0) {
  stop("No Sample/Subgenome set in the gene-copy CSV matches the phase-log columns.")
}

expected_logs <- paste0(prefix, "_locus_", seq_along(loci), "_phase.log")
if (any(!file.exists(expected_logs))) {
  stop("The number/order of locus columns in the gene-copy CSV does not match the phase logs.")
}

# =============================================================================
# 5. FUNCTION THAT ADDS THE EXACT-TIP HEATMAP TO THE TREE
# =============================================================================
# Modified from ggtree::gheatmap. It draws only assignments observed in the
# joint-MAP state. Outlined white boxes retain the alignment of BLANK or
# unreported cells. Only focal constrained tips receive heatmap boxes.
add_exact_tip_heatmap <- function(
    p, probability_data, label_data, offset = 0, width = 1,
    palette = "cividis", colnames_position = "top", font_size = 6,
    legend_title = "Marginal posterior probability", width_reference = NULL) {
  colnames_position <- match.arg(colnames_position, c("bottom", "top"))

  tree_data <- p$data
  if (is.null(width_reference)) {
    width_reference <- diff(range(tree_data$x, na.rm = TRUE))
  }
  heatmap_width <- width * width_reference / ncol(probability_data)
  cell_width <- heatmap_width * 0.87

  node_coordinates <- intersect(
    unlist(tree_data %>% filter(is.na(.data$x)) %>% select(.data$parent, .data$node)),
    unlist(tree_data %>% filter(!is.na(.data$x)) %>% select(.data$parent, .data$node))
  )
  selected_internal_labels <- tree_data %>%
    filter(.data$node %in% node_coordinates) %>%
    select(.data$label) %>%
    unlist()
  selected_labels <- intersect(selected_internal_labels, rownames(probability_data))
  keep <- tree_data$isTip | tree_data$label %in% selected_labels
  display_tree_data <- tree_data[keep, ]
  start <- max(display_tree_data$x, na.rm = TRUE) + offset

  row_order <- order(display_tree_data$y)
  row_order <- row_order[!is.na(display_tree_data$y[row_order])]
  labels_in_tree_order <- display_tree_data$label[row_order]

  probabilities <- as.data.frame(probability_data)
  labels <- as.data.frame(label_data, stringsAsFactors = FALSE)
  probabilities <- probabilities[
    match(labels_in_tree_order, rownames(probabilities)), , drop = FALSE
  ]
  labels <- labels[match(labels_in_tree_order, rownames(labels)), , drop = FALSE]

  probabilities$y <- sort(display_tree_data$y)
  labels$y <- sort(display_tree_data$y)
  probabilities$lab <- labels_in_tree_order
  labels$lab <- labels_in_tree_order

  long_probabilities <- tidyr::pivot_longer(
    probabilities, cols = -c("lab", "y"),
    names_to = "variable", values_to = "value"
  )
  long_labels <- tidyr::pivot_longer(
    labels, cols = -c("lab", "y"),
    names_to = "variable", values_to = "display_label"
  )
  long_probabilities$value <- as.numeric(long_probabilities$value)
  long_probabilities$display_label <- as.character(long_labels$display_label)
  long_probabilities$display_label[
    is.na(long_probabilities$display_label) |
      !nzchar(long_probabilities$display_label) |
      grepl("BLANK", long_probabilities$display_label, ignore.case = TRUE)
  ] <- NA_character_
  long_probabilities$value[is.na(long_probabilities$display_label)] <- NA_real_
  long_probabilities$variable <- factor(
    long_probabilities$variable, levels = colnames(probability_data)
  )
  long_probabilities$x <- start +
    as.numeric(long_probabilities$variable) * heatmap_width
  long_probabilities[[".panel"]] <- factor("Tree")

  column_mapping <- unique(data.frame(
    from = long_probabilities$variable,
    to = long_probabilities$x,
    stringsAsFactors = FALSE
  ))
  column_mapping[[".panel"]] <- factor("Tree")
  column_mapping$y <- if (colnames_position == "bottom") {
    0
  } else {
    max(tree_data$y) + 1
  }

  tiles <- long_probabilities[
    !is.na(long_probabilities$value) &
      !is.na(long_probabilities$display_label), , drop = FALSE
  ]

  # Draw the complete locus grid only for tips that belong to the focal
  # analysis. Reference-tree tips stay untouched. BLANK/unreported cells are
  # therefore visible as empty white boxes rather than disappearing.
  grid_tiles <- long_probabilities[
    long_probabilities$lab %in% rownames(probability_data), , drop = FALSE
  ]

  p2 <- p + geom_tile(
    data = grid_tiles,
    aes(x = .data$x, y = .data$y),
    height = 0.81, width = cell_width,
    fill = "white", color = NA, inherit.aes = FALSE
  )
  p2 <- p2 + geom_tile(
    data = tiles,
    aes(x = .data$x, y = .data$y, fill = .data$value),
    height = 0.81, width = cell_width, color = NA, inherit.aes = FALSE
  )
  p2 <- p2 + geom_tile(
    data = grid_tiles,
    aes(x = .data$x, y = .data$y),
    height = 0.81, width = cell_width,
    fill = NA, color = blank_box_border_color,
    linewidth = blank_box_border_width, inherit.aes = FALSE
  )
  p2 <- p2 + geom_text(
    data = tiles[tiles$value >= 0.58, , drop = FALSE],
    aes(x = .data$x, y = .data$y, label = .data$display_label),
    size = 3.7, fontface = "bold", color = "black", inherit.aes = FALSE
  )
  p2 <- p2 + geom_text(
    data = tiles[tiles$value < 0.58, , drop = FALSE],
    aes(x = .data$x, y = .data$y, label = .data$display_label),
    size = 3.7, fontface = "bold", color = "white", inherit.aes = FALSE
  )

  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    exact_colorbar <- guide_colorbar(
      direction = "vertical",
      theme = theme(
        legend.title.position = "left",
        legend.title = element_text(
          angle = 90, face = "bold", size = 10.5,
          hjust = 0.5, vjust = 0.5
        ),
        legend.key.height = grid::unit(2.4, "in"),
        legend.key.width = grid::unit(0.4125, "in"),
        legend.text = element_text(size = 9, face = "bold")
      )
    )
  } else {
    exact_colorbar <- guide_colorbar(
      direction = "vertical", title.position = "left",
      title.theme = element_text(
        angle = 90, face = "bold", size = 10.5,
        hjust = 0.5, vjust = 0.5
      ),
      title.hjust = 0.5, title.vjust = 0.5,
      barheight = grid::unit(2.4, "in"),
      barwidth = grid::unit(0.4125, "in"),
      label.theme = element_text(size = 9, face = "bold")
    )
  }

  p2 <- p2 + scale_fill_viridis_c(
    option = palette, direction = 1, limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1), na.value = "white",
    name = legend_title, guide = exact_colorbar
  )
  p2 <- p2 + geom_text(
    data = column_mapping,
    aes(x = .data$to, y = .data$y, label = .data$from),
    size = font_size, fontface = "plain", inherit.aes = FALSE
  )

  legend_theme <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.03, 0.95),
      legend.justification = c(0, 1)
    )
  } else {
    theme(
      legend.position = c(0.03, 0.95),
      legend.justification = c(0, 1)
    )
  }
  p2 <- p2 + legend_theme + theme(
    legend.background = element_rect(
      fill = grDevices::adjustcolor("white", alpha.f = 0.95),
      color = "grey70", linewidth = 0.25
    ),
    legend.margin = margin(6, 9, 6, 7.5)
  )

  attr(p2, "mapping") <- column_mapping
  p2
}

# =============================================================================
# 6. CALCULATE JOINT-MAP ASSIGNMENTS AND EXACT-TIP MARGINAL PP
# =============================================================================
# For each locus and sample, the script first identifies the most frequently
# sampled complete phase arrangement after burn-in. For every displayed copy,
# exact-tip marginal PP is then the fraction of retained samples in which that
# copy occupies that exact constrained tip.

all_tips <- unique(unlist(samples, use.names = FALSE))
exact_tip_probabilities <- as.data.frame(
  matrix(
    NA_real_, nrow = length(all_tips), ncol = length(loci),
    dimnames = list(all_tips, loci)
  )
)
joint_map_phase_results <- as.data.frame(
  matrix(
    "", nrow = length(all_tips), ncol = length(loci),
    dimnames = list(all_tips, loci)
  ),
  stringsAsFactors = FALSE
)
exact_tip_copy_confidence <- data.frame()

for (locus_index in seq_along(loci)) {
  phase_log <- paste0(prefix, "_locus_", locus_index, "_phase.log")
  all_draws <- read.delim(
    phase_log, stringsAsFactors = FALSE, check.names = FALSE
  )
  first_post_burnin <- floor(nrow(all_draws) * burnin) + 1
  draws_after_burnin <- all_draws[
    first_post_burnin:nrow(all_draws), , drop = FALSE
  ]

  for (sample in names(samples)) {
    sample_tips <- samples[[sample]]
    draws <- draws_after_burnin[, sample_tips, drop = FALSE]

    # The displayed copy arrangement is the most frequently sampled complete
    # assignment. Ties are resolved deterministically by table ordering.
    state_key <- apply(draws, 1, paste, collapse = "\u001f")
    state_frequency <- sort(table(state_key), decreasing = TRUE)
    map_row <- which(state_key == names(state_frequency)[1])[1]
    map_phases <- as.character(draws[map_row, sample_tips])
    names(map_phases) <- sample_tips
    joint_map_phase_results[sample_tips, loci[locus_index]] <- map_phases

    for (tip in sample_tips) {
      copy <- map_phases[tip]
      is_reported <- !is.na(copy) && nzchar(copy) &&
        !grepl("BLANK", copy, ignore.case = TRUE)
      if (!is_reported) {
        exact_tip_probabilities[tip, loci[locus_index]] <- NA_real_
        next
      }

      exact_pp <- mean(draws[[tip]] == copy, na.rm = TRUE)
      exact_tip_probabilities[tip, loci[locus_index]] <- exact_pp
      exact_tip_copy_confidence <- rbind(
        exact_tip_copy_confidence,
        data.frame(
          sample = sample,
          locus = loci[locus_index],
          tip = tip,
          copy = copy,
          exact_tip_marginal_pp = exact_pp,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# =============================================================================
# 7. BUILD AND SAVE THE FIGURE
# =============================================================================
# Branch lengths are multiplied only in the in-memory tree object. The original
# tree file and every posterior probability remain unchanged.

for (palette_name in palette_names) {
  tree <- treeio::read.beast(tree_file)
  if (is.null(tree@phylo$edge.length)) {
    stop("The input tree has no branch lengths to scale.")
  }
  original_tree_span <- max(ape::node.depth.edgelength(tree@phylo))
  tree@phylo$edge.length <- tree@phylo$edge.length * branch_length_scale
  plot_height <- max(10, ape::Ntip(tree@phylo) * 0.33)
  p <- ggtree(tree)
  p <- p + geom_tiplab(
    size = 4.55, fontface = "bold", align = TRUE,
    linesize = 0.3, offset = 0.0005
  )
  p <- add_exact_tip_heatmap(
    p, exact_tip_probabilities, joint_map_phase_results,
    offset = 0.020, width = 0.95, palette = palette_name,
    colnames_position = "top", font_size = 6,
    legend_title = "Marginal posterior probability",
    width_reference = original_tree_span
  )
  p <- p + theme(
    legend.text = element_text(size = 9, face = "bold"),
    legend.title = element_text(
      angle = 90, size = 10.5, face = "bold",
      hjust = 0.5, vjust = 0.5
    ),
    axis.text = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )
  p <- p + ggtree::geom_nodelab(
    geom = "text", nudge_x = -0.00099, nudge_y = 0.25,
    size = 3.9, fontface = "bold",
    ggplot2::aes(label = round(as.numeric(posterior), 2))
  )

  file_stem <- paste0(
    "tree_plot_exact_tip_marginal_pp_", palette_name, "_", output_label
  )
  ggsave(
    file.path(figure_dir, paste0(file_stem, ".pdf")),
    plot = p, width = 16, height = plot_height, units = "in"
  )
  ggsave(
    file.path(figure_dir, paste0(file_stem, ".png")),
    plot = p, width = 16, height = plot_height, units = "in", dpi = 300
  )
}

# =============================================================================
# 8. SAVE MACHINE-READABLE RESULTS
# =============================================================================
# The confidence table contains one row per displayed copy-tip assignment. The
# joint-MAP table is the phase arrangement used for the text inside the cells.

write.csv(
  exact_tip_copy_confidence,
  file = file.path(
    figure_dir, paste0("exact_tip_copy_confidence_summary_", output_label, ".csv")
  ),
  row.names = FALSE
)
write.csv(
  joint_map_phase_results,
  file = file.path(
    figure_dir, paste0("joint_map_phase_results_", output_label, ".csv")
  )
)

message("Completed exact-tip marginal PP plot for: ", paste(names(samples), collapse = ", "))
message("Uniform branch-length scale: ", branch_length_scale)
message("Detected output directory: ", output_dir)
message("Results written to: ", normalizePath(figure_dir, winslash = "/"))
