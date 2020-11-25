![](https://static.inaturalist.org/photos/64272938/original.jpg?1584917392)
# Usage Notes
This repository contains scripts and data required to replicate the analyses.

## Raw Data
Raw, unprocessed sequencing data is provided in `*.fastq.gz` format un the folder `1_raw_data`. The processed data (012NA genotype files) are in the folder `PAWR_WIWR_012NA_files`. These data include the genotype files with the suffix `.012NA` (contain genotypes for each sample at each position), position files with the suffix `.pos` (listing chromosome name and position for each SNP), and individual files with the suffix `.indv` (listing the sample order of each individual in the dataset).

## Sample metadata
Sample metadata are provided in the folder `Sample_information`.

## Scripts
Instructions for running each step of the analysis is provided in text files in markdown format:  
1) Instructions for processing the raw data is provided in `1_process_sequences.md`
2) Instructions for performing PCA on the data is provided in `2_PCA_plotting.md` and associated data are in the folder `2_PCA_data`
3) Instructions for the STRUCTURE analysis is provided in `3_STRUCTURE_analysis.md` and associated raw data are in the folder `3_STRUCTURE_data`
4) Instructions for the phylogenetic network are provided in `4_Phylogenetic_network.md` and associated data are in the folder `4_SplitsTree`
5) Instructions for performing the genome scans (Fst, pi_within, and pi_between) are in `5_GenomeScans.md`
