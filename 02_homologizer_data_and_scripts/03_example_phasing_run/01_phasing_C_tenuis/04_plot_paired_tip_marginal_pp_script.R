# =============================================================================
# Paired-tip marginal posterior-probability plot for Homologizer
# =============================================================================
#
# PURPOSE
# -------
# This script summarizes a paired-tip constraint analysis at the biological
# subgenome/homeolog-pair level. If a copy moves between B1 and B2 across MCMC
# samples, the plot reports the probability that it occupies either member of
# the B pair rather than treating B1 and B2 as different assignments.
#
# A single displayed copy produces one cell spanning the two paired rows. If
# two copies are assigned to the same pair, the cell is divided vertically and
# each copy retains its own pair-marginal posterior probability (PP).
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
# Paired tips must end in 1 and 2, for example A1/A2 or B1/B2.

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

# Destination for the PDF, PNG, and two CSV summaries. "auto" creates a
# paired_tip_plot_results folder inside analysis_dir. You may instead supply a
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
  figure_dir <- file.path(analysis_dir, "paired_tip_plot_results")
} else {
  figure_dir <- path.expand(figure_dir)
  figure_dir_is_absolute <- grepl(
    "^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", figure_dir
  )
  if (!figure_dir_is_absolute) {
    figure_dir <- file.path(analysis_dir, figure_dir)
  }
}

allowed_palettes <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")
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

# Auto-detect the Homologizer output directory when requested.
if (identical(tolower(output_folder), "auto")) {
  # Accept either the analysis parent directory or the Homologizer output
  # directory itself, which makes batch use and ad-hoc reruns equally simple.
  candidate_dirs <- unique(c(
    analysis_dir,
    list.dirs(analysis_dir, recursive = FALSE, full.names = TRUE)
  ))
  candidate_dirs <- candidate_dirs[vapply(candidate_dirs, function(path) {
    length(list.files(path, pattern = "^homologizer_locus_[0-9]+_phase[.]log$")) > 0
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

# Locate the gene-copy map by its required Sample and Subgenome columns. Search
# the analysis directory and a few ancestors so the input can be either a
# temp_outputs directory or output1/output2 itself.
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
  header <- tryCatch(names(read.csv(path, nrows = 1, check.names = FALSE)), error = function(e) character())
  all(c("Sample", "Subgenome") %in% header)
}, logical(1))]
if (length(valid_csv) == 0) {
  stop(
    "Could not find a gene-copy mapping CSV containing Sample and Subgenome ",
    "columns (normally subgenomes_names.csv)."
  )
}
genecopyFn <- valid_csv[1]

# Prefer map1 when present, otherwise use map, then any unique map tree.
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
prefix <- file.path(output_dir, "homologizer")

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

genecopymap = read.csv(genecopyFn, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
samples = split(as.character(genecopymap$Subgenome), genecopymap$Sample)

# names of the loci in the log file
loci = names(genecopymap)[3:length(genecopymap)]

# Retain only samples whose complete paired-tip set is present in the logs.
first_log <- paste0(prefix, "_locus_1_phase.log")
if (!file.exists(first_log)) {
  stop("Missing expected log file: ", first_log)
}
log_columns <- names(read.delim(first_log, nrows = 1, check.names = FALSE))
samples <- samples[vapply(samples, function(tips) all(tips %in% log_columns), logical(1))]
if (length(samples) == 0) {
  stop("No Sample/Subgenome set in the gene-copy CSV matches the phase-log columns.")
}

expected_logs <- paste0(prefix, "_locus_", seq_along(loci), "_phase.log")
if (any(!file.exists(expected_logs))) {
  stop("The number/order of locus columns in the gene-copy CSV does not match the phase logs.")
}

# =============================================================================
# 5. FUNCTION THAT ADDS THE PAIRED-TIP HEATMAP TO THE TREE
# =============================================================================
# This function is modified from ggtree::gheatmap. It constructs pair-spanning
# cells, vertical divisions for two-copy pairs, and outlined white boxes for
# BLANK/unreported assignments. Only focal constrained tips receive boxes.

add_paired_tip_heatmap = function (p, data, data_labels, offset = 0, width = 1,
                        palette = "cividis",
                        merge_single_copy_pairs = FALSE,
                        show_row_summary = TRUE,
                        color = "white", colnames = TRUE, colnames_position = "bottom",
                        colnames_angle = 0, colnames_level = NULL, colnames_offset_x = 0,
                        colnames_offset_y = 0, font.size = 4, family = "", hjust = 0.5,
                        legend_title = "Marginal posterior probability",
                        width_reference = NULL)
{
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  if (is.null(width_reference)) {
    width_reference <- p$data$x %>% range(na.rm = TRUE) %>% diff
  }
  width <- width * width_reference/ncol(data)
  # Keep a visible gutter between locus columns. The paired-copy divider is
  # drawn separately as a much finer line.
  cell_width <- width * 0.87
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  df <- p$data
  nodeCo <- intersect(df %>% filter(is.na(x)) %>% select(.data$parent,
                                                         .data$node) %>% unlist(), df %>% filter(!is.na(x)) %>%
                        select(.data$parent, .data$node) %>% unlist())
  labCo <- df %>% filter(.data$node %in% nodeCo) %>% select(.data$label) %>%
    unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- df$label %in% selCo
  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm = TRUE) + offset
  dd <- as.data.frame(data)
  dd2 <- as.data.frame(data_labels)
  i <- order(df$y)
  i <- i[!is.na(df$y[i])]
  lab <- df$label[i]
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  dd2 <- dd2[match(lab, rownames(dd2)), , drop = FALSE]
  dd$y <- sort(df$y)
  dd2$y <- sort(df$y)
  dd$lab <- lab
  dd2$lab <- lab
  dd <- gather(dd, variable, value, -c(lab, y))
  dd2 <- gather(dd2, variable, value, -c(lab, y))
  dd$value <- as.numeric(dd$value)
  dd2$value[dd2$value == ""] <- NA
  dd2$probability <- dd$value
  dd$display_label <- as.character(dd2$value)
  dd$tile_height <- 0.81
  dd$tile_width <- cell_width
  dd$divider_x <- NA_real_
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  }
  else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from = dd$variable, to = V2)
  mapping <- unique(mapping)
  dd$x <- V2
  dd2$x <- V2
  dd$width <- width
  dd2$width <- width
  dd[[".panel"]] <- factor("Tree")
  dd2[[".panel"]] <- factor("Tree")

  # Create a complete white background grid only for focal constrained tips.
  # In the paired/homeolog view, one empty box spans both paired rows, matching
  # the geometry used by a reported one-copy homeolog assignment.
  grid_source <- dd[dd$lab %in% rownames(data), , drop = FALSE]
  if (merge_single_copy_pairs) {
    grid_source$pair_group <- homeolog_group(grid_source$lab)
    grid_key <- interaction(
      grid_source$variable, grid_source$pair_group, drop = TRUE
    )
    grid_parts <- lapply(split(grid_source, grid_key), function(group_data) {
      if (nrow(group_data) == 2 && diff(range(group_data$y)) <= 1.01) {
        blank_cell <- group_data[1, , drop = FALSE]
        blank_cell$y <- mean(group_data$y)
        blank_cell$tile_height <- diff(range(group_data$y)) + 0.81
        return(blank_cell)
      }
      group_data
    })
    grid_data <- do.call(rbind, grid_parts)
  } else {
    grid_data <- grid_source
  }

  # Draw only reported assignments. In the homeolog-level plot, a single
  # observed copy assigned to an A1/A2, B1/B2, C1/C2, or D1/D2 pair is drawn as
  # one tile spanning both rows. If two observed copies occupy the pair, the
  # same pair-spanning cell is divided left/right with a vertical boundary;
  # each half carries that copy's own pair-collapsed marginal PP.
  if (merge_single_copy_pairs) {
    dd$pair_group <- homeolog_group(dd$lab)
    split_key <- interaction(dd$variable, dd$pair_group, drop = TRUE)
    tile_parts <- lapply(split(dd, split_key), function(group_data) {
      reported <- which(
        !is.na(group_data$value) & !is.na(group_data$display_label)
      )
      if (length(reported) == 0) {
        return(group_data[FALSE, , drop = FALSE])
      }
      if (
        nrow(group_data) == 2 && length(reported) == 1 &&
          diff(range(group_data$y)) <= 1.01
      ) {
        merged <- group_data[reported, , drop = FALSE]
        merged$y <- mean(group_data$y)
        merged$tile_height <- diff(range(group_data$y)) + 0.81
        return(merged)
      }
      if (
        nrow(group_data) == 2 && length(reported) == 2 &&
          diff(range(group_data$y)) <= 1.01
      ) {
        vertical_split <- group_data[reported, , drop = FALSE]
        vertical_split <- vertical_split[order(vertical_split$y), , drop = FALSE]
        vertical_split$y <- mean(group_data$y)
        vertical_split$tile_height <- diff(range(group_data$y)) + 0.81
        vertical_split$tile_width <- cell_width / 2
        vertical_split$divider_x <- unique(group_data$x)[1]
        vertical_split$x <- unique(group_data$x)[1] +
          c(-cell_width / 4, cell_width / 4)
        return(vertical_split)
      }
      group_data[reported, , drop = FALSE]
    })
    tile_data <- do.call(rbind, tile_parts)
  }
  else {
    tile_data <- dd[!is.na(dd$value) & !is.na(dd$display_label), , drop = FALSE]
  }

  p2 <- p + geom_tile(
    data = grid_data,
    aes(x, y, height = tile_height, width = tile_width),
    fill = "white", color = NA, inherit.aes = FALSE
  )

  if (is.null(color)) {
    p2 <- p2 + geom_tile(
      data = tile_data,
      aes(x, y, fill = value, height = tile_height, width = tile_width),
      inherit.aes = FALSE
    )
  }
  else {
    p2 <- p2 + geom_tile(
      data = tile_data,
      aes(x, y, fill = value, height = tile_height, width = tile_width),
      color = NA, inherit.aes = FALSE
    )
    tile_data$is_vertical_split <- tile_data$tile_width < cell_width * 0.75
    divider_data <- unique(tile_data[
      tile_data$is_vertical_split,
      c("divider_x", "y", "tile_height")
    ])
    if (nrow(divider_data) > 0) {
      p2 <- p2 + geom_segment(
        data = divider_data,
        aes(
          x = divider_x, xend = divider_x,
          y = y - tile_height / 2, yend = y + tile_height / 2
        ),
        color = "white", linewidth = 0.2, inherit.aes = FALSE
      )
    }
    # Switch label color with luminance so low-confidence dark tiles remain
    # legible. Use slightly smaller type in the two narrow vertical halves.
    p2 <- p2 + geom_text(
      data = tile_data[tile_data$value >= 0.58 & !tile_data$is_vertical_split, ],
      aes(x, y, label = display_label), size = 3.7, fontface = "bold",
      color = "black", inherit.aes = FALSE
    )
    p2 <- p2 + geom_text(
      data = tile_data[tile_data$value < 0.58 & !tile_data$is_vertical_split, ],
      aes(x, y, label = display_label), size = 3.7, fontface = "bold",
      color = "white", inherit.aes = FALSE
    )
    p2 <- p2 + geom_text(
      data = tile_data[tile_data$value >= 0.58 & tile_data$is_vertical_split, ],
      aes(x, y, label = display_label), size = 2.62, fontface = "bold",
      color = "black", lineheight = 0.9, inherit.aes = FALSE
    )
    p2 <- p2 + geom_text(
      data = tile_data[tile_data$value < 0.58 & tile_data$is_vertical_split, ],
      aes(x, y, label = display_label), size = 2.62, fontface = "bold",
      color = "white", lineheight = 0.9, inherit.aes = FALSE
    )

    if (show_row_summary) {
      dd3 = data.frame()
      start_x = max(dd$x)
      height = max(dd$y)
      margin = 0.006
      for (y in unique(tile_data$y)) {
        pp = mean(tile_data[tile_data$y == y, 'value'], na.rm=TRUE)
        dd4 = data.frame(pp = pp, x = pp/200 + margin + start_x, y=y)
        dd3 = rbind(dd3, dd4)
      }
      p2 <- p2 + geom_segment(aes(x=start_x+margin, xend=1/200 + start_x + margin, y=0.2, yend=0.2), size=0.5, inherit.aes = FALSE)
      p2 <- p2 + geom_segment(aes(x=1/200+start_x+margin, xend=1/200+start_x+margin, y=0.5, yend=height), color='grey85', linetype='dotted', size=0.35, inherit.aes = FALSE)
      p2 <- p2 + geom_segment(aes(x=start_x+margin, xend=start_x+margin, y=0.5, yend=height), color='grey85', linetype='dotted', size=0.35, inherit.aes = FALSE)
      p2 <- p2 + geom_point(data = dd3, aes(x, y, color=pp), size=1.25, inherit.aes = FALSE, show.legend=FALSE)
      p2 <- p2 + geom_text(label='0.0', x=start_x+margin, y=-0.2, size=1.25, color='grey50')
      p2 <- p2 + geom_text(label='1.0', x=start_x+margin+1/200, y=-0.2, size=1.25, color='grey50')
    }
  }
  # Draw outlines last so both filled and empty cells have the same crisp outer
  # boundary. Blank cells remain white and contain no text.
  p2 <- p2 + geom_tile(
    data = grid_data,
    aes(x, y, height = tile_height, width = tile_width),
    fill = NA, color = blank_box_border_color,
    linewidth = blank_box_border_width, inherit.aes = FALSE
  )
  if (methods::is(dd$value, "numeric")) {
    # A monotonic, color-vision-deficiency-friendly scale makes 0.1, 0.5, and
    # 1.0 visually ordered. Turbo is retained only for direct comparison.
    # Build the legend guide according to the installed ggplot2 version. This
    # keeps the title vertical and immediately left of the bar in both the
    # modern guide-theme API and the older title.position/title.theme API.
    if (utils::packageVersion("ggplot2") >= "3.5.0") {
      homeolog_colorbar <- guide_colorbar(
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
      homeolog_colorbar <- guide_colorbar(
        direction = "vertical",
        title.position = "left",
        title.theme = element_text(
          angle = 90, face = "bold", size = 10.5,
          hjust = 0.5, vjust = 0.5
        ),
        title.hjust = 0.5,
        title.vjust = 0.5,
        barheight = grid::unit(2.4, "in"),
        barwidth = grid::unit(0.4125, "in"),
        label.theme = element_text(size = 9, face = "bold")
      )
    }
    p2 <- p2 + scale_fill_viridis_c(
      option = palette, direction = 1, limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      na.value = "white", name = legend_title,
      guide = homeolog_colorbar
    )
    p2 <- p2 + scale_color_viridis_c(
      option = palette, direction = 1, limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      na.value = "grey80", name = legend_title
    )
  }
  else {
    p2 <- p2 + scale_fill_discrete(na.value = NA, name = legend_title)
  }
  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    }
    else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    p2 <- p2 + geom_text(data = mapping, aes(x = to, y = y,
                                             label = from), size = font.size, family = family,
                         inherit.aes = FALSE, angle = colnames_angle, nudge_x = colnames_offset_x,
                         nudge_y = colnames_offset_y, hjust = hjust)
  }
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
  if (!colnames) {
    p2 <- p2 + scale_y_continuous(expand = c(0, 0))
  }
  attr(p2, "mapping") <- mapping
  return(p2)
}



# =============================================================================
# 6. IDENTIFY THE CONSTRAINED TIP PAIRS
# =============================================================================
# Infer valid exchangeable pairs from the actual tips. A trailing 1 or 2 is
# collapsed only when both members exist, so unrelated names ending in a digit
# are not accidentally merged.
make_pair_groups <- function(tips) {
  candidate <- sub("[12]$", "", tips)
  suffix <- ifelse(grepl("[12]$", tips), sub(".*([12])$", "\\1", tips), NA_character_)
  groups <- tips
  for (base in unique(candidate)) {
    members <- which(candidate == base & !is.na(suffix))
    if (length(members) == 2 && setequal(suffix[members], c("1", "2"))) {
      groups[members] <- base
    }
  }
  groups
}

pair_lookup <- do.call(c, unname(lapply(samples, function(tips) {
  setNames(make_pair_groups(tips), tips)
})))
homeolog_group <- function(x) {
  matched <- match(x, names(pair_lookup))
  result <- x
  result[!is.na(matched)] <- unname(pair_lookup[matched[!is.na(matched)]])
  result
}

pair_definition <- do.call(rbind, lapply(names(samples), function(sample) {
  tips <- samples[[sample]]
  data.frame(
    sample = sample,
    tip = tips,
    pair_group = homeolog_group(tips),
    stringsAsFactors = FALSE
  )
}))
paired_counts <- table(pair_definition$pair_group)
if (!any(paired_counts == 2)) {
  stop("No valid terminal 1/2 tip pairs were detected.")
}

all_tips <- unique(unlist(samples, use.names = FALSE))
homeolog_prob_results <- as.data.frame(
  matrix(NA_real_, nrow = length(all_tips), ncol = length(loci),
         dimnames = list(all_tips, loci))
)
joint_map_phase_results <- as.data.frame(
  matrix("", nrow = length(all_tips), ncol = length(loci),
         dimnames = list(all_tips, loci)),
  stringsAsFactors = FALSE
)
homeolog_phase_results <- joint_map_phase_results
copy_confidence_summary <- data.frame()

# =============================================================================
# 7. SUMMARIZE THE POSTERIOR SAMPLES
# =============================================================================
# For each locus and focal accession, the script first finds the joint MAP
# assignment. It then calculates, for every displayed copy, the fraction of
# post-burn-in samples in which that copy occurs on either member of its pair.

for (i in seq_along(loci)) {
  f_in <- paste0(prefix, "_locus_", i, "_phase.log")
  d_all <- read.delim(f_in, stringsAsFactors = FALSE, check.names = FALSE)
  first_post_burnin <- floor(nrow(d_all) * burnin) + 1
  d_all <- d_all[first_post_burnin:nrow(d_all), , drop = FALSE]

  for (sample in names(samples)) {
    sample_tips <- samples[[sample]]
    sample_groups <- homeolog_group(sample_tips)
    d <- d_all[, sample_tips, drop = FALSE]

    # Serialize each sampled joint state instead of building a potentially huge
    # Cartesian-product table of every possible phase combination.
    state_key <- apply(d, 1, paste, collapse = "\u001f")
    state_frequency <- sort(table(state_key), decreasing = TRUE)
    map_row <- which(state_key == names(state_frequency)[1])[1]
    map_phases <- as.character(d[map_row, sample_tips])
    names(map_phases) <- sample_tips
    joint_map_phase_results[sample_tips, loci[i]] <- map_phases

    observed <- !is.na(map_phases) & nzchar(map_phases) &
      !grepl("BLANK", map_phases, ignore.case = TRUE)

    for (tip_index in seq_along(sample_tips)) {
      tip <- sample_tips[tip_index]
      phase <- map_phases[tip]
      exact_pp <- mean(d[[tip]] == phase, na.rm = TRUE)

      if (!observed[tip_index]) {
        homeolog_prob_results[tip, loci[i]] <- NA_real_
        homeolog_phase_results[tip, loci[i]] <- ""
        next
      }

      group_tips <- sample_tips[sample_groups == sample_groups[tip_index]]
      collapsed_pp <- mean(apply(
        d[, group_tips, drop = FALSE], 1, function(state) phase %in% state
      ), na.rm = TRUE)
      homeolog_prob_results[tip, loci[i]] <- collapsed_pp
      homeolog_phase_results[tip, loci[i]] <- phase
      copy_confidence_summary <- rbind(
        copy_confidence_summary,
        data.frame(
          sample = sample,
          locus = loci[i],
          homeolog = sample_groups[tip_index],
          copy = phase,
          displayed_tip = tip,
          exact_tip_marginal_pp = exact_pp,
          homeolog_assignment_pp = collapsed_pp,
          stringsAsFactors = FALSE
        )
      )
    }

  }
}

# =============================================================================
# 8. BUILD AND SAVE THE PAIRED-TIP FIGURE
# =============================================================================
# branch_length_scale changes only this in-memory copy of the tree. The original
# tree span is retained as the reference for heatmap widths and offsets, so those
# elements remain visually stable when the branches are made more compact.
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
  p <- add_paired_tip_heatmap(
    p, homeolog_prob_results, homeolog_phase_results,
    offset = 0.020, palette = palette_name,
    merge_single_copy_pairs = TRUE,
    show_row_summary = FALSE,
    colnames_position = "top", font.size = 6, width = 0.95,
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
    geom = "text",
    nudge_x = -0.00099,
    nudge_y = 0.25,
    size = 3.9,
    fontface = "bold",
    ggplot2::aes(label = round(as.numeric(posterior), 2))
  )

  file_stem <- paste0(
    "tree_plot_paired_tip_marginal_pp_", palette_name, "_", output_label
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
# 9. SAVE THE TWO REUSABLE SUMMARY TABLES
# =============================================================================
write.csv(
  joint_map_phase_results,
  file = file.path(figure_dir, paste0("joint_map_phase_results_", output_label, ".csv"))
)
write.csv(
  copy_confidence_summary,
  file = file.path(figure_dir, paste0("copy_confidence_summary_", output_label, ".csv")),
  row.names = FALSE
)

message("Completed paired-tip marginal-PP plot for: ", paste(names(samples), collapse = ", "))
message("Uniform branch-length scale: ", branch_length_scale)
message("Detected output directory: ", output_dir)
message("Results written to: ", normalizePath(figure_dir, winslash = "/"))
