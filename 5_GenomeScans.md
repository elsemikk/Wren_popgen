This file details the pipeline used to perform genome scans of Fst, pi_within, and pi_between for the Pacific and Winter Wren dataset.

**Input**: `*.012NA`, `*.indv`, and `*.pos` files containing the genotype data for all samples ,produced in the first step. These files are provided in the folder `2_PCA_DATA` in this repository. It also requires metadata files which are provided in the same folder. This R code also uses functions that are provided in the R script `genomics_R_functions.R`.

**Output**: estimates of Fst, pi_within, and pi_between, for windows across the genome.

