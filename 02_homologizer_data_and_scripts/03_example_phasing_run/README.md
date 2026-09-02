# Worked homologizer analyses

These directories contain the locus alignments, RevBayes inputs, and derived
figures/tables for representative analyses in the manuscript. Together they
show the progression from a straightforward paired-tip analysis to the
diagnosis of allelic ambiguity, the phasing of a complex taxon, and a
fixed-topology backbone check.

## Example index

| Directory | Analysis | What it demonstrates |
|---|---|---|
| [`01_phasing_C_tenuis/`](01_phasing_C_tenuis/) | Paired-tip phasing of *Cystopteris tenuis* 16972 | A relatively simple tetraploid example used during iterative backbone construction |
| [`02_Comparison_of_C_deqinensis_phasing/`](02_Comparison_of_C_deqinensis_phasing/) | Standard and paired-tip analyses of *C. deqinensis* 12509 | Why allelic copies can make the standard setup ambiguous and how paired tips resolve the interpretation (Fig. 4) |
| [`03_phasing_C_laurentiana/`](03_phasing_C_laurentiana/) | Paired-tip phasing of *C. laurentiana* 7886 | Application of the backbone to a more complex polyploid (Fig. 5) |
| [`04_backbone_rephasing_fixed_tree/`](04_backbone_rephasing_fixed_tree/) | Rephasing backbone polyploids on a fixed tree | Consistency check underlying the phase summary in Fig. 3 |
| [`05_scripts_used_to_phase_other_accessions/`](05_scripts_used_to_phase_other_accessions/) | Accession-specific RevBayes input files | Full record of the iterative analyses used to assemble the phased backbone |

## *Cystopteris deqinensis* comparison

The comparison is divided into three matched analyses:

1. `01_standard_2tips/`: one tip per expected subgenome;
2. `02_standard_3tips/`: an additional unconstrained tip accommodates a second
   allelic copy; and
3. `03_paired_tips/`: two constrained sister tips represent each expected
   biological subgenome.

The exact-tip plot for the paired analysis may show intermediate values because
a copy can exchange between the two members of a constrained pair. The
paired-tip plot integrates over that within-pair exchange and reports support
for assignment to the biological subgenome represented by the pair.

## Contents of each example

- `data/`: FASTA alignments used by RevBayes.
- `homologizer_scripts/`: model, initial phase, phase moves, and any topology or
  paired-tip constraints.
- `subgenomes_name.csv` or `subgenomes_names.csv`: focal tip map used by the R
  plotting scripts.
- `exact_tip_plot_results/`: exact-tip PDF/PNG and CSV summaries.
- `paired_tip_plot_results/`: paired-tip PDF/PNG and CSV summaries when the
  analysis used paired constraints.

The plotting scripts located directly inside examples are synchronized copies
of the canonical scripts in
[`../02_homologizer_paired_tip_phasing_pipeline/`](../02_homologizer_paired_tip_phasing_pipeline/).
For adaptation and future maintenance, use the canonical copies.

## Rerunning an example

1. Enter the example's `homologizer_scripts/` directory.
2. Create the sibling output directory specified by `output_file` in the main
   RevBayes script, usually `../output/`.
3. Start RevBayes and source the example's main `.Rev` script.
4. Assess convergence and summarize the posterior.
5. Open the appropriate canonical R plotting script, set `analysis_dir` to the
   example directory and `output_folder` to the newly generated output folder,
   then run the entire R script.

Precomputed figures and summary CSVs are included for inspection. Raw MCMC
traces are not retained in every example, so they must be regenerated before a
plot can be reproduced from phase logs.
