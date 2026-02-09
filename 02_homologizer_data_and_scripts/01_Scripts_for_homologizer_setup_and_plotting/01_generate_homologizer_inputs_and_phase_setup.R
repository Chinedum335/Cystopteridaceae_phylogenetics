# ============================================================
# Script: 01_generate_homologizer_inputs.R
# Project: Cystopteridaceae Review – MUL-tree Analyses
# Author: Chinedum Anajemba
# Description:
#   Generates key input files for homologizer and RevBayes:
#     - subgenome_names.csv
#     - paired_tip_constraint.rev
#     - InitialPhase.rev
#     - PhaseMoves.rev
#
# Input:
#   Generate_phase.csv
#
# Output:
#   Files written to working directory or /revbayes/
# ============================================================


# Load libraries
library(readr)
library(dplyr)
library(stringr)

# Read the input CSV
book2 <- read.csv("Generate_phase.csv")
df <- read.csv("Generate_phase.csv")

# ---- PhaseMaker Input Creation ----
generate_rows <- function(sample, sample_id, no_subgenome, app, gap, ibr, pgi, trngr) {
  n_upper <- ceiling(no_subgenome / 2)
  n_lower <- floor(no_subgenome / 2)
  
  # Instead of plain A/a, use A1, A2
  upper_suffixes <- if (n_upper > 0) paste0(LETTERS[1:n_upper], "1") else character()
  lower_suffixes <- if (n_lower > 0) paste0(LETTERS[1:n_lower], "2") else character()
  
  suffixes <- c(upper_suffixes, lower_suffixes)
  subgenomes <- paste0(sample, "_", suffixes)
  
  create_column_data <- function(count, no_subgenome, sample_id) {
    if (count == 0) {
      return(paste0(sample_id, "_BLANK", seq_len(no_subgenome)))
    } else {
      copies <- paste0(sample_id, "_c", 1:count)
      blanks <- paste0(sample_id, "_BLANK", seq_len(no_subgenome - count))
      return(c(copies, blanks)[1:no_subgenome])
    }
  }
  
  data.frame(
    Sample = rep(sample, no_subgenome),
    Subgenome = subgenomes,
    APP = create_column_data(app, no_subgenome, sample_id),
    GAP = create_column_data(gap, no_subgenome, sample_id),
    IBR = create_column_data(ibr, no_subgenome, sample_id),
    PGI = create_column_data(pgi, no_subgenome, sample_id),
    TRNGR = create_column_data(trngr, no_subgenome, sample_id),
    stringsAsFactors = FALSE
  )
}

book1 <- data.frame()
for (i in 1:nrow(book2)) {
  row <- book2[i, ]
  book1 <- bind_rows(book1, generate_rows(
    sample = row$Sample,
    sample_id = row$sample_ID,
    no_subgenome = row$No_subgenome,
    app = row$APP,
    gap = row$GAP,
    ibr = row$IBR,
    pgi = row$PGI,
    trngr = row$TRNGR
  ))
}


cystopteridaceae_genomes <- book1 %>%
  mutate(APP = "", GAP = "", IBR = "", PGI = "", TRNGR = "")
write.csv(cystopteridaceae_genomes, "subgenomes_names.csv", row.names = FALSE)

# ---- Clade Constraint File ----
unique_samples <- unique(book1$Sample)
rev_lines <- c()
all_clade_names <- c()

for (sample in unique_samples) {
  sample_id <- unique(book2$sample_ID[book2$Sample == sample])
  subgenomes <- book1 %>% filter(Sample == sample) %>% pull(Subgenome)
  
  # Loop through LETTERS and look for A1/A2, B1/B2, etc.
  for (letter in LETTERS) {
    sub1 <- paste0(sample, "_", letter, "1")
    sub2 <- paste0(sample, "_", letter, "2")
    if (sub1 %in% subgenomes & sub2 %in% subgenomes) {
      clade_name <- paste0("clade_", sample_id, "_", letter)
      clade_def <- sprintf('%s = clade("%s", "%s")', clade_name, sub1, sub2)
      rev_lines <- c(rev_lines, clade_def)
      all_clade_names <- c(all_clade_names, clade_name)
    }
  }
}

rev_lines <- c(
  rev_lines,
  "",
  "# combine the constraints",
  sprintf("constraints = [%s]", paste(all_clade_names, collapse = ", "))
)

writeLines(rev_lines, "paired_tip_constraint.rev")

# ---- InitialPhase.rev and PhaseMoves.rev ----
genecopymap <- book1
numLoci <- ncol(genecopymap) - 2

missingtaxa.text <- "data[%d].addMissingTaxa(\"%s\")\n"
initial.phase.text <- "data[%d].setHomeologPhase(\"%s\", \"%s\")\n"
mv.text <- "moves[++mvi] = mvHomeologPhase(ctmc[%d], \"%s\", \"%s\", weight=3)\n"

missingtaxa.commands <- c()
initial.phase.commands <- c()
my.index <- 1
for (gene in 3:ncol(genecopymap)) {
  for (sample in 1:nrow(genecopymap)) {
    val <- as.character(genecopymap[sample, gene])
    if (grepl("BLANK", val)) {
      missingtaxa.commands <- c(missingtaxa.commands, sprintf(missingtaxa.text, gene-2, val))
    }
    initial.phase.commands[my.index] <- sprintf(initial.phase.text, gene-2, val, genecopymap[sample, "Subgenome"])
    my.index <- my.index + 1
  }
}

cat(missingtaxa.commands, file = "InitialPhase.rev")
cat(initial.phase.commands, file = "InitialPhase.rev", append = TRUE)

subgenomes.by.sample <- split(genecopymap$Subgenome, genecopymap$Sample)
make.phase.comb <- function(subgenome.samples, geneNum) {
  out_text <- ""
  for (i in 1:(length(subgenome.samples) - 1)) {
    for (j in (i + 1):length(subgenome.samples)) {
      out_text <- paste0(out_text, sprintf(mv.text, geneNum, subgenome.samples[i], subgenome.samples[j]))
    }
  }
  return(out_text)
}

mv.commands <- c()
for (gene in 1:numLoci) {
  for (sample in 1:length(subgenomes.by.sample)) {
    mv.commands <- c(mv.commands, make.phase.comb(subgenomes.by.sample[[sample]], gene))
  }
}
cat(mv.commands, file = "PhaseMoves.rev")
