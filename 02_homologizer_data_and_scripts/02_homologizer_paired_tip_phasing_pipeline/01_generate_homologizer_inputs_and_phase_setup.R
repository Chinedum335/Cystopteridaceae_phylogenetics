# =============================================================================
# Generate Homologizer phase files for the paired-tip constraint approach
# =============================================================================
#
# PURPOSE
# -------
# Starting from Generate_phase.csv, this script creates the four files needed
# to describe focal polyploid samples and their exchangeable paired tips:
#
#   subgenomes_names.csv      mapping table used by the plotting scripts
#   paired_tip_constraint.rev sister-tip constraints used by RevBayes
#   InitialPhase.rev          starting copy-to-tip assignments
#   PhaseMoves.rev            Homologizer phase-swap proposals
#
# HOW TO RUN IN RSTUDIO
# ---------------------
# 1. Put this script and Generate_phase.csv in the same working directory.
# 2. Edit the USER SETTINGS section if necessary.
# 3. Select the entire script and run it, or click Source.
#
# This implementation uses base R only.

# =============================================================================
# 1. USER SETTINGS
# =============================================================================

phase_config_file <- "Generate_phase.csv"
output_dir <- "."
locus_columns <- c("APP", "GAP", "IBR", "PGI", "TRNGR")

# =============================================================================
# 2. READ AND VALIDATE THE CONFIGURATION
# =============================================================================

if (!file.exists(phase_config_file)) {
  stop("Cannot find phase configuration file: ", phase_config_file)
}

phase_config <- read.csv(
  phase_config_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("Sample", "sample_ID", "No_subgenome", locus_columns)
missing_columns <- setdiff(required_columns, names(phase_config))
if (length(missing_columns) > 0) {
  stop(
    "Generate_phase.csv is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
if (nrow(phase_config) == 0) {
  stop("Generate_phase.csv contains no samples.")
}
if (any(!nzchar(trimws(phase_config$Sample))) ||
    any(!nzchar(trimws(as.character(phase_config$sample_ID))))) {
  stop("Sample and sample_ID must be supplied for every row.")
}

numeric_columns <- c("No_subgenome", locus_columns)
for (column in numeric_columns) {
  phase_config[[column]] <- suppressWarnings(as.numeric(phase_config[[column]]))
  if (any(!is.finite(phase_config[[column]])) ||
      any(phase_config[[column]] != floor(phase_config[[column]]))) {
    stop(column, " must contain whole numbers only.")
  }
}

if (any(phase_config$No_subgenome < 2) ||
    any(phase_config$No_subgenome %% 2 != 0)) {
  stop(
    "No_subgenome must be a positive even number because every biological ",
    "subgenome is represented by two constrained tips."
  )
}
for (locus in locus_columns) {
  if (any(phase_config[[locus]] < 0) ||
      any(phase_config[[locus]] > phase_config$No_subgenome)) {
    stop(
      locus, " copy counts must be between 0 and No_subgenome for each sample."
    )
  }
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir)
}
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

# =============================================================================
# 3. BUILD THE INITIAL COPY-TO-TIP TABLE
# =============================================================================

make_copy_vector <- function(copy_count, tip_count, sample_id) {
  copies <- if (copy_count > 0) {
    paste0(sample_id, "_c", seq_len(copy_count))
  } else {
    character(0)
  }
  blank_count <- tip_count - copy_count
  blanks <- if (blank_count > 0) {
    paste0(sample_id, "_BLANK", seq_len(blank_count))
  } else {
    character(0)
  }
  c(copies, blanks)
}

make_sample_rows <- function(config_row) {
  tip_count <- config_row$No_subgenome
  pair_count <- tip_count / 2
  pair_letters <- LETTERS[seq_len(pair_count)]

  # Ordering A1, B1, ..., A2, B2, ... preserves the historical setup while
  # terminal 1/2 labels identify the two exchangeable tips in each pair.
  suffixes <- c(paste0(pair_letters, "1"), paste0(pair_letters, "2"))
  result <- data.frame(
    Sample = rep(config_row$Sample, tip_count),
    Subgenome = paste0(config_row$Sample, "_", suffixes),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (locus in locus_columns) {
    result[[locus]] <- make_copy_vector(
      config_row[[locus]], tip_count, as.character(config_row$sample_ID)
    )
  }
  result
}

phase_rows <- do.call(rbind, lapply(seq_len(nrow(phase_config)), function(i) {
  make_sample_rows(phase_config[i, , drop = FALSE])
}))
rownames(phase_rows) <- NULL

# The plotting map needs tip names but not the arbitrary initial assignments.
plotting_map <- phase_rows
plotting_map[locus_columns] <- ""
write.csv(
  plotting_map,
  file.path(output_dir, "subgenomes_names.csv"),
  row.names = FALSE
)

# =============================================================================
# 4. WRITE THE PAIRED-TIP CLADE CONSTRAINTS
# =============================================================================

constraint_lines <- character(0)
constraint_names <- character(0)
for (i in seq_len(nrow(phase_config))) {
  sample_name <- phase_config$Sample[i]
  sample_id <- as.character(phase_config$sample_ID[i])
  pair_count <- phase_config$No_subgenome[i] / 2
  for (pair_letter in LETTERS[seq_len(pair_count)]) {
    constraint_name <- paste0("clade_", sample_id, "_", pair_letter)
    tip_1 <- paste0(sample_name, "_", pair_letter, "1")
    tip_2 <- paste0(sample_name, "_", pair_letter, "2")
    constraint_lines <- c(
      constraint_lines,
      sprintf('%s = clade("%s", "%s")', constraint_name, tip_1, tip_2)
    )
    constraint_names <- c(constraint_names, constraint_name)
  }
}
constraint_lines <- c(
  constraint_lines,
  "",
  "# Combine all paired-tip constraints.",
  sprintf("constraints = [%s]", paste(constraint_names, collapse = ", "))
)
writeLines(
  constraint_lines,
  file.path(output_dir, "paired_tip_constraint.rev")
)

# =============================================================================
# 5. WRITE INITIALPHASE.REV
# =============================================================================

missing_taxon_lines <- character(0)
initial_phase_lines <- character(0)
for (locus_index in seq_along(locus_columns)) {
  locus <- locus_columns[locus_index]
  for (row_index in seq_len(nrow(phase_rows))) {
    copy_name <- phase_rows[[locus]][row_index]
    tip_name <- phase_rows$Subgenome[row_index]
    if (grepl("BLANK", copy_name, fixed = TRUE)) {
      missing_taxon_lines <- c(
        missing_taxon_lines,
        sprintf('data[%d].addMissingTaxa("%s")', locus_index, copy_name)
      )
    }
    initial_phase_lines <- c(
      initial_phase_lines,
      sprintf(
        'data[%d].setHomeologPhase("%s", "%s")',
        locus_index, copy_name, tip_name
      )
    )
  }
}
writeLines(
  c(missing_taxon_lines, initial_phase_lines),
  file.path(output_dir, "InitialPhase.rev")
)

# =============================================================================
# 6. WRITE PHASEMOVES.REV
# =============================================================================

tips_by_sample <- split(phase_rows$Subgenome, phase_rows$Sample)
phase_move_lines <- character(0)
for (locus_index in seq_along(locus_columns)) {
  for (sample_tips in tips_by_sample) {
    tip_pairs <- utils::combn(sample_tips, 2)
    for (pair_index in seq_len(ncol(tip_pairs))) {
      phase_move_lines <- c(
        phase_move_lines,
        sprintf(
          'moves[++mvi] = mvHomeologPhase(ctmc[%d], "%s", "%s", weight=3)',
          locus_index, tip_pairs[1, pair_index], tip_pairs[2, pair_index]
        )
      )
    }
  }
}
writeLines(phase_move_lines, file.path(output_dir, "PhaseMoves.rev"))

message("Created paired-tip setup files in: ", output_dir)
message("Samples: ", paste(phase_config$Sample, collapse = ", "))
message("Loci: ", paste(locus_columns, collapse = ", "))
