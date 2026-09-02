# homologizer data, paired-tip workflow, and worked analyses

This directory contains the phased Cystopteridaceae backbone sequences, the
general paired-tip homologizer workflow, and worked analyses from the
manuscript.

## Directory guide

| Directory | Contents | Use it when... |
|---|---|---|
| [`01_Backbone_phased_sequences/`](01_Backbone_phased_sequences/) | Five locus-specific FASTA files with established subgenome assignments encoded in sequence names | You want a backbone reference for phasing additional Cystopteridaceae accessions |
| [`02_homologizer_paired_tip_phasing_pipeline/`](02_homologizer_paired_tip_phasing_pipeline/) | General R and RevBayes scripts for generating, running, and plotting paired-tip analyses | You want to adapt the approach to a new accession or another organismal system |
| [`03_example_phasing_run/`](03_example_phasing_run/) | Worked standard-tip, paired-tip, complex-polyploid, and fixed-tree analyses | You want concrete input/output examples or the analyses associated with manuscript figures |

## Backbone sequences

The files in `01_Backbone_phased_sequences/data/` are the final phased backbone
alignments for APP, GAP, IBR, PGI, and TRNGR. Sequence names encode the inferred
subgenome assignment. For example, if copy `16972_c1` was assigned to the A
subgenome of *Cystopteris tenuis* 16972, the fixed backbone sequence carries the
corresponding `_A` tip identity.

When reusing this backbone:

1. add the phased backbone sequences to the corresponding locus alignments for
   the new focal accession;
2. retain the backbone tip names exactly;
3. treat backbone sequences as fixed/reference tips; and
4. assign `mvHomeologPhase` proposals only to copies from the new focal
   accession(s).

## Iterative construction of the backbone

The backbone was assembled using a diploids-first, iterative phasing strategy.
Polyploids with better-understood ancestry were phased first. Their inferred
subgenome assignments were then fixed by renaming their sequences, and those
fixed sequences were added as reference tips in subsequent analyses. This was
repeated until the backbone contained representatives of the major inferred
subgenome histories.

The scripts retained in
`03_example_phasing_run/05_scripts_used_to_phase_other_accessions/` document the
accession-specific analyses used during that process.

## Standard versus paired-tip analyses

A standard homologizer analysis uses one tip for each expected subgenome. The
paired-tip workflow uses two constrained sister tips for each expected
subgenome. The latter can accommodate two allelic sequences belonging to one
subgenome without forcing them to represent distinct homeologous histories.

The two plotting scripts report different marginal probabilities:

- exact-tip plot: probability that a copy occupies one named tip;
- paired-tip plot: probability that a copy occupies either member of the
  constrained sister pair representing a subgenome.

See the [pipeline README](02_homologizer_paired_tip_phasing_pipeline/README.md)
for the full workflow and interpretation.
