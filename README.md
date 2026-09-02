# Cystopteridaceae phylogenetics and paired-tip homologizer workflow

This repository contains the data, RevBayes scripts, R utilities, worked
examples, and derived plotting outputs associated with the manuscript
*Cystopteridaceae in focus: a model clade for studying polyploid systematics and
macroevolution*.

The repository has two linked purposes:

1. document the PURC and homologizer analyses used to build and phase the
   Cystopteridaceae backbone MUL-tree; and
2. provide an adaptable paired-tip workflow for distinguishing allelic copies
   within a subgenome from homeologous copies belonging to different
   subgenomes.

## Start here

| Goal | Go to |
|---|---|
| Inspect the PURC configuration, barcodes, and locus maps | [`01_purc/`](01_purc/) |
| Use the phased Cystopteridaceae backbone sequences | [`01_Backbone_phased_sequences/`](02_homologizer_data_and_scripts/01_Backbone_phased_sequences/) |
| Apply the paired-tip approach to a new dataset | [`02_homologizer_paired_tip_phasing_pipeline/`](02_homologizer_data_and_scripts/02_homologizer_paired_tip_phasing_pipeline/) |
| Reproduce or inspect the manuscript examples | [`03_example_phasing_run/`](02_homologizer_data_and_scripts/03_example_phasing_run/) |
| Plot exact-tip marginal probabilities | [`05_plot_exact_tip_marginal_pp_script.R`](02_homologizer_data_and_scripts/02_homologizer_paired_tip_phasing_pipeline/05_plot_exact_tip_marginal_pp_script.R) |
| Plot paired-tip/subgenome marginal probabilities | [`04_plot_paired_tip_marginal_pp_script.R`](02_homologizer_data_and_scripts/02_homologizer_paired_tip_phasing_pipeline/04_plot_paired_tip_marginal_pp_script.R) |

For a new paired-tip analysis, begin with the [pipeline guide](02_homologizer_data_and_scripts/02_homologizer_paired_tip_phasing_pipeline/README.md). For the biological and file-level organization of the homologizer material, see the [homologizer data guide](02_homologizer_data_and_scripts/README.md).

## Repository layout

```text
.
|-- 01_purc/
|   |-- README.md
|   |-- PURC_config_file.txt
|   |-- Barcode_seq.fasta
|   |-- CY25_*map.txt
|   `-- Reference_seq.fa
|-- 02_homologizer_data_and_scripts/
|   |-- README.md
|   |-- 01_Backbone_phased_sequences/
|   |   `-- data/
|   |-- 02_homologizer_paired_tip_phasing_pipeline/
|   |   |-- README.md
|   |   |-- 01_generate_homologizer_inputs_and_phase_setup.R
|   |   |-- 02_homologizer_main_script.Rev
|   |   |-- 03_homologizer_main_script_no_clade_constraint.Rev
|   |   |-- 04_plot_paired_tip_marginal_pp_script.R
|   |   |-- 05_plot_exact_tip_marginal_pp_script.R
|   |   `-- Generate_phase.csv
|   `-- 03_example_phasing_run/
|       |-- README.md
|       |-- 01_phasing_C_tenuis/
|       |-- 02_Comparison_of_C_deqinensis_phasing/
|       |-- 03_phasing_C_laurentiana/
|       |-- 04_backbone_rephasing_fixed_tree/
|       `-- 05_scripts_used_to_phase_other_accessions/
|-- .gitignore
`-- README.md
```

The numbered directory names preserve the order of the analytical workflow.
Within a worked example, the usual layout is:

```text
example/
|-- data/                  locus-specific FASTA alignments
|-- homologizer_scripts/   RevBayes model and generated phase files
|-- exact_tip_plot_results/
|-- paired_tip_plot_results/   paired-tip examples only
`-- subgenomes_names.csv   focal-sample/tip mapping used for plotting
```

## Conceptual overview

### Standard homologizer phasing

[homologizer](https://doi.org/10.1111/2041-210X.14072) jointly estimates a
phylogeny and the assignment of recovered gene copies to the tips representing
the evolutionary histories of a polyploid. In a standard analysis, one tip is
usually supplied for each expected subgenome.

When two sequences recovered from a locus are allelic variants from the same
subgenome, a one-tip-per-subgenome model can force those sequences onto
different subgenome tips. Adding an extra tip can diagnose this situation, but
doing so requires recognizing the ambiguity and rerunning the analysis.

### Paired-tip constraint approach

The paired-tip approach represents every expected biological subgenome with two
constrained sister tips. For example, subgenome B is represented by B1 and B2.
This gives two allelic copies from the same subgenome somewhere to coexist while
the pair remains a single constrained lineage in the tree.

For a tetraploid expected to contain two biological subgenomes, the input
therefore contains four phase tips: A1, A2, B1, and B2. The constraint file
requires A1/A2 and B1/B2 to be sisters. homologizer can still exchange copies
among all phase tips during MCMC.

## Interpreting the two posterior-probability plots

The two plotting scripts answer related but different questions.

| Plot | Quantity displayed | Recommended use |
|---|---|---|
| Exact-tip marginal PP | Probability that a copy is assigned to one exact tip, such as B1 | Standard single-tip analyses; diagnosis of how probability is divided within a paired-tip analysis |
| Paired-tip marginal PP | Probability that a copy is assigned to either member of a pair, such as B1 or B2 | Primary presentation of paired-tip/subgenome assignments |

In a paired-tip run, B1 and B2 are exchangeable. A copy that is
always assigned to subgenome B may spend about half of the posterior samples on
B1 and half on B2. The exact-tip marginal PP would then be approximately 0.5,
even though support for assignment to the B pair is approximately 1.0. The
paired-tip plotting script sums over that within-pair uncertainty and reveals
the biologically relevant subgenome-level support.

If only one copy is assigned to a pair at a locus, the paired-tip plot uses one
box spanning the two tip rows. If two copies are assigned to the same pair, the
box is divided vertically so that each copy and its own pair-level marginal PP
remain visible. White outlined boxes indicate loci for which no assignment is
displayed.

## Software requirements

### RevBayes analysis

- [RevBayes](https://revbayes.github.io/) with homologizer functions available,
  including `mvHomeologPhase` and `mnHomeologPhase`
- A command-line shell or the RevBayes interactive console

### R setup and plotting

- R
- CRAN packages: `ape`, `dplyr`, `ggplot2`, `magrittr`, and `tidyr`
- Bioconductor packages: `ggtree` and `treeio`

Install the plotting dependencies from R with:

```r
install.packages(c("ape", "dplyr", "ggplot2", "magrittr", "tidyr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ggtree", "treeio"))
```

The phase-file generator uses base R and has no additional package dependency.

## Quick start: adapt the paired-tip workflow

### 1. Create an analysis directory

```text
my_analysis/
|-- data/
|-- homologizer_scripts/
`-- output/
```

Copy the five locus alignments into `data/`. Copy
`01_generate_homologizer_inputs_and_phase_setup.R`,
`02_homologizer_main_script.Rev`, and `Generate_phase.csv` into
`homologizer_scripts/`.

The supplied RevBayes script assumes it is run from `homologizer_scripts/`,
reads alignments from `../data/`, and writes files with the prefix
`../output/homologizer`.

### 2. Name recovered copies consistently

Copies from each focal accession must use the accession identifier followed by
`_c1`, `_c2`, and so on. For example:

```text
7886_c1
7886_c2
7886_c3
```

Reference/backbone sequences should already have fixed tip names and should not
receive phase moves.

### 3. Describe the focal samples

Edit `Generate_phase.csv`. It has one row per focal accession:

| Column | Meaning |
|---|---|
| `Sample` | Base taxon/accession name used to construct tip labels |
| `sample_ID` | Prefix used by the recovered copies in the FASTA files |
| `No_subgenome` | Total number of phase tips; this must be twice the number of expected biological subgenomes |
| `APP`–`TRNGR` | Number of copies recovered at each locus |

The included example describes *Cystopteris tenuis* 16972 with two biological
subgenomes, represented by four constrained phase tips.

### 4. Generate the phase and constraint files

In RStudio, set the working directory to `homologizer_scripts/`, open
`01_generate_homologizer_inputs_and_phase_setup.R`, select the complete script,
and run it. Alternatively:

```bash
cd my_analysis/homologizer_scripts
Rscript 01_generate_homologizer_inputs_and_phase_setup.R
```

The script validates the configuration and creates:

- `InitialPhase.rev`
- `PhaseMoves.rev`
- `paired_tip_constraint.rev`
- `subgenomes_names.csv`

### 5. Run RevBayes

From `homologizer_scripts/`, open RevBayes and source the paired-tip model:

```text
source("02_homologizer_main_script.Rev")
```

The template currently performs 1,000 tuning/burn-in generations followed by
15,000 main generations and summarizes the posterior trees after discarding 30%
of the tree trace. These settings are visible near the end of the Rev script and
should be evaluated for convergence in every new analysis.

### 6. Plot the results

Open the desired plotting script in RStudio and edit only its **USER SETTINGS**
section. For example:

```r
analysis_dir <- "C:/path/to/my_analysis"
output_folder <- "output"
figure_dir <- "auto"
palette_names <- c("cividis")
burnin <- 0.20
branch_length_scale <- 1.00
```

Select the entire script and run it. `figure_dir <- "auto"` creates a results
directory inside `analysis_dir`; a full destination path can be supplied
instead. `branch_length_scale = 1.00` preserves the original tree branch
lengths, while values below 1 make the plotted tree more compact without
changing the input tree.

## Plotting outputs

Each plotting script writes:

- a PDF figure;
- a 300-dpi PNG figure;
- the joint MAP phase assignment table; and
- a copy-level posterior-probability summary.

The paired-tip copy summary includes both the exact-tip marginal PP and the
paired-tip/subgenome marginal PP. Both figure legends are titled **Marginal
posterior probability**; the assignment being marginalized differs as described
above.

## Worked examples and manuscript figures

The [worked-example guide](02_homologizer_data_and_scripts/03_example_phasing_run/README.md)
maps each directory to its scientific purpose. In particular:

- *C. deqinensis* compares standard two-tip, standard three-tip, and paired-tip
  analyses and underlies the methodological comparison in Fig. 4 of the manuscript.
- *C. laurentiana* demonstrates use of the paired-tip backbone to phase a more
  complex taxon and underlies Fig. 5.
- the fixed-tree backbone rephasing analysis corresponds to the phase summary
  in Fig. 3.

Precomputed plot PDFs, PNGs, and CSV summaries are retained with the examples.
The large raw MCMC traces are not included in every example directory; rerun the
provided RevBayes scripts to regenerate them before rerunning a plot locally.

## Sequence datasets

The repository uses five loci:

- nuclear: APP (ApPEFP-C), GAP (gapCpSh), IBR (IBR3), and PGI (pgiC)
- plastid: TRNGR

The phased backbone FASTA files encode inferred subgenome assignments directly
in their sequence names. When those sequences are reused as a reference for a
new analysis, treat them as fixed tips rather than assigning new phasing moves.

## Reproducibility notes

- Run each RevBayes example from its `homologizer_scripts/` directory because
  paths in the model files are relative to that location.
- Assess convergence across independent chains (run about 6 independent chains); the supplied generation counts
  are analysis settings, not a universal guarantee of adequate sampling.
- Keep the row order and locus order consistent among the configuration file,
  alignments, generated phase files, and plotting map.
- R session files (`.Rhistory` and `.RData`) are ignored because they are not
  reproducible analysis inputs.

## Related resources and citation

- Freyman, W. A., Johnson, M. G., and Rothfels, C. J. (2023).
  [homologizer: Phylogenetic phasing of gene copies into polyploid subgenomes](https://doi.org/10.1111/2041-210X.14072).
  *Methods in Ecology and Evolution*.
- [Original homologizer example repository](https://github.com/wf8/homologizer)
- [RevBayes documentation](https://revbayes.github.io/)
- Rothfels et al. (2017). [Next-generation polyploid phylogenetics: rapid resolution of hybrid polyploid complexes using PacBio single-molecule sequencing](https://doi.org/10.1111/nph.14111).
- Schafran et al. (2023). [PURC Provides Improved Sequence Inference for Polyploid Phylogenetics and Other Manifestations of the Multiple-Copy Problem] (https://link.springer.com/protocol/10.1007/978-1-0716-2561-3_10).
- Project data archive: [Dryad](https://doi.org/10.5061/dryad.brv15dvqc)

When citing this repository, please cite the accompanying Cystopteridaceae
article and the archived dataset once their final bibliographic details are
available.

## Data availability

The repository is maintained at
[github.com/Chinedum335/Cystopteridaceae_phylogenetics](https://github.com/Chinedum335/Cystopteridaceae_phylogenetics).
Final sequence alignments, scripts, and voucher metadata are also archived in
[Dryad](https://doi.org/10.5061/dryad.brv15dvqc).

## License

No explicit software or data license is currently included. Until the authors
select one, reuse is governed by ordinary copyright and the terms of the Dryad
record. Adding a repository license before publication is recommended.
