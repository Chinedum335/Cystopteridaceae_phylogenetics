# PURC input files

This directory contains the configuration and reference files used to process
the long-read amplicon data with PURC (Pipeline for Untangling Reticulate
Complexes).

## Contents

| File | Purpose |
|---|---|
| `PURC_config_file.txt` | Analysis configuration used for PURC processing |
| `Barcode_seq.fasta` | Pacific Biosciences barcode sequences used for demultiplexing |
| `CY25_APPmap.txt` | Sample/barcode map for APP |
| `CY25_GAPmap.txt` | Sample/barcode map for GAP |
| `CY25_IBRmap.txt` | Sample/barcode map for IBR |
| `CY25_PGImap.txt` | Sample/barcode map for PGI |
| `CY25_TRNGRmap.txt` | Sample/barcode map for TRNGR |
| `Reference_seq.fa` | Cystopteridaceae reference sequences used during processing |

These are study-specific inputs, not a complete installation of PURC. Install
PURC and its dependencies separately, then provide these files in the locations
expected by the PURC configuration. Before adapting the configuration, inspect
every path and parameter because local directory paths and computing resources
will differ among systems.

For the method and software background, see Rothfels et al. (2017) and Schafran et al. (2023),
[Next-generation polyploid phylogenetics](https://doi.org/10.1111/nph.14111).
[PURC Provides Improved Sequence Inference for Polyploid Phylogenetics](https://link.springer.com/protocol/10.1007/978-1-0716-2561-3_10).
