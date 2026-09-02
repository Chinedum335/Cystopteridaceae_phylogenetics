# Paired-tip homologizer pipeline

This directory is the canonical, reusable implementation of the paired-tip
constraint approach. The scripts can be copied into a new analysis and adapted
by editing the clearly marked input and settings sections.

## Files

| File | Role |
|---|---|
| `Generate_phase.csv` | Describes focal samples, total phase tips, and recovered copy counts by locus |
| `01_generate_homologizer_inputs_and_phase_setup.R` | Validates `Generate_phase.csv` and creates phase, move, constraint, and plotting-map files |
| `02_homologizer_main_script.Rev` | Runs joint phylogeny and phase inference with paired-tip topological constraints |
| `03_homologizer_main_script_no_clade_constraint.Rev` | Runs the standard homologizer model without paired-tip constraints |
| `04_plot_paired_tip_marginal_pp_script.R` | Recommended plot for paired-tip results; summarizes support across both members of each subgenome pair |
| `05_plot_exact_tip_marginal_pp_script.R` | Exact-tip plot for standard analyses or for diagnosing within-pair exchangeability |

## Expected analysis layout

The supplied RevBayes scripts use relative paths. A portable analysis should be
organized as follows:

```text
my_analysis/
|-- data/
|   |-- APP_homologizer.fa
|   |-- GAP_homologizer.fa
|   |-- IBR_homologizer.fa
|   |-- PGI_homologizer.fa
|   `-- TRNGR_homologizer.fa
|-- homologizer_scripts/
|   |-- Generate_phase.csv
|   |-- 01_generate_homologizer_inputs_and_phase_setup.R
|   |-- 02_homologizer_main_script.Rev
|   |-- InitialPhase.rev
|   |-- PhaseMoves.rev
|   |-- paired_tip_constraint.rev
|   `-- subgenomes_names.csv
`-- output/
```

Run RevBayes from `homologizer_scripts/`. The paired model reads the five FASTA
files from `../data/` and writes the analysis under `../output/`.

## Step 1: prepare alignments

Each locus alignment must contain the fixed reference/backbone tips and the
unphased copies from each focal sample. Focal copies must share the identifier
listed in `sample_ID` and use `_c1`, `_c2`, etc. For example:

```text
12509_c1
12509_c2
```

The locus columns and their order must agree among the FASTA files,
`Generate_phase.csv`, the R generator, and the RevBayes `alignments` vector.

## Step 2: edit `Generate_phase.csv`

Required columns:

| Column | Description |
|---|---|
| `Sample` | Base name used for the new phylogenetic tips |
| `sample_ID` | Prefix of the recovered copies in the FASTA files |
| `No_subgenome` | Total phase tips, not the biological subgenome count; it must be a positive even number |
| Locus columns | Number of copies recovered for each locus |

For two biological subgenomes, use `No_subgenome = 4`; the generated tips are
A1, B1, A2, and B2, and the constraint pairs are A1/A2 and B1/B2. For three
biological subgenomes, use `No_subgenome = 6`.

## Step 3: generate phase inputs

In RStudio, set the working directory to `homologizer_scripts/`, open
`01_generate_homologizer_inputs_and_phase_setup.R`, select the entire script,
and run it. From a shell:

```bash
cd my_analysis/homologizer_scripts
Rscript 01_generate_homologizer_inputs_and_phase_setup.R
```

The generator creates:

- `InitialPhase.rev`: an arbitrary valid starting assignment, including
  placeholder taxa for unrecovered copies;
- `PhaseMoves.rev`: pairwise homologizer phase-swap moves among focal phase
  tips;
- `paired_tip_constraint.rev`: sister constraints for the two tips representing
  each biological subgenome; and
- `subgenomes_names.csv`: focal sample/tip map read by the plotting scripts.

The starting assignment does not assert biological certainty. Phase assignments
are sampled by the MCMC.

## Step 4: run the paired-tip model

Create the `output/` directory, start RevBayes from `homologizer_scripts/`, and
run:

```text
source("02_homologizer_main_script.Rev")
```

The template jointly samples topology, branch lengths, GTR+Gamma substitution
parameters, among-locus rates, and homeolog phase. It writes a model log, tree
trace, locus-specific phase logs, MAP tree, and MCC tree.

The template uses 1,000 tuning/burn-in generations followed by 15,000 main
generations. The MAP/MCC tree summary discards 30% of the tree trace. Treat
these as documented starting settings: convergence and effective sample sizes
must be assessed for each dataset and across independent runs.

## Step 5: plot posterior assignments

Open either plotting script in RStudio and edit its **USER SETTINGS** section:

```r
analysis_dir <- "C:/path/to/my_analysis"
output_folder <- "output"
figure_dir <- "auto"
palette_names <- c("cividis")
burnin <- 0.20
branch_length_scale <- 1.00
```

Then select the entire script and run it. The scripts locate the MAP tree,
locus-specific phase logs, and a mapping CSV containing `Sample` and
`Subgenome`. They produce PDF/PNG figures and two CSV summaries.

### Which plot should be reported?

Use `04_plot_paired_tip_marginal_pp_script.R` as the primary summary of a
paired-tip analysis. It reports the posterior probability that a copy is
assigned to the biological subgenome represented by either member of a tip
pair. Use `05_plot_exact_tip_marginal_pp_script.R` for standard single-tip
analyses or to show how paired-tip support is divided between exchangeable
members.

Both legends read **Marginal posterior probability**. In the exact-tip plot the
event is assignment to one tip; in the paired-tip plot the event is assignment
to either tip in the pair.

### Display conventions

- Cell text is the copy assigned in the joint maximum a posteriori phase state.
- Cell color is the relevant marginal posterior probability for that copy.
- One copy assigned to a pair produces one box spanning the paired rows.
- Two copies assigned to the same pair divide that box vertically.
- White outlined cells mean no assignment is displayed.
- `branch_length_scale = 1.00` preserves the input branch lengths. Values below
  1 compact the tree for display only.

## Adapting beyond Cystopteridaceae

The paired-tip setup is not taxon-specific. To adapt it:

1. replace the reference and focal sequences in each locus alignment;
2. preserve unique copy identifiers using a consistent `_cN` convention;
3. change `locus_columns` in the generator and the `alignments` vector in both
   RevBayes models if the number or names of loci differ;
4. ensure the mapping CSV has `Sample`, `Subgenome`, and matching locus columns;
5. retain terminal `1` and `2` on paired-tip names because the paired plotting
   script uses those suffixes to identify exchangeable pairs; and
6. adjust MCMC length, priors, and proposal weights to the new dataset and
   verify convergence.

## Troubleshooting

- **Multiple output folders found:** set `output_folder` explicitly rather than
  using `"auto"`.
- **Cannot find mapping CSV:** keep `subgenomes_names.csv` in the analysis
  directory or a nearby parent directory.
- **Cannot write the figure:** use `figure_dir <- "auto"` or provide a shorter,
  writable absolute path; close an existing PDF if it is locked.
- **Messages about an invalid edge matrix or deprecated `.data`:** these can be
  emitted by compatible but differently versioned `ggtree`, `treeio`, and
  `tidytree` packages. They are warnings, not the output-writing error.
- **Paired script reports no valid pairs:** confirm that both members exist and
  their tip names end in `1` and `2`.
