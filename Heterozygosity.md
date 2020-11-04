


```R
# set directory where files stored
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")

# Load functions
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/genomics_R_functions.R") #differs from old code at line 205

# choose the chromosomes to analyze in this run
chromosomes.to.analyze <- c("autoinfSNP")  

Analysis_set <- 1  # 1: all samples, only SNPs;    2: all samples, with invariant sites

if (Analysis_set == 1) {
  # choose path and filename for the 012NA files
  base.name <- "PAWR_WIWR_012NA_files/PAWR_WIWR.genotypes.SNPs_only"
  filename.text.middle <- ".max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 30 
  filename.text.end <-".MQ20.lowHet.tab"
  tag.name <- ".all_chr."   # choose a tag name for this analysis
  # indicate name of metadata file, a text file with these column headings:
  # ID  location  group Fst_group plot_order
  metadata.file <- "PAWRWIWR_Metadata_autoinfSNP.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 75  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR", "Hybrid")
  group.colors <- c("blue", "red", "purple")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR")
  group.colors.WC84_Fst <- c("purple")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR", "Hybrid")
  group.colors.pi <- c("blue", "red", "purple")
  
} else if (Analysis_set == 2) {
  # choose path and filename for the 012NA files
  base.name <- "PAWR_WIWR_012NA_files/PAWR_WIWR.genotypes.allSites"
  filename.text.middle <- ".infoSites.max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 30 
  filename.text.end <-".MQ20.lowHet.tab"
  tag.name <- ".all_samples_with_invariant."   # choose a tag name for this analysis
  # indicate name of metadata file, a text file with these column headings:
  # ID  location  group Fst_group plot_order
  metadata.file <- "PAWRWIWR_Metadata_autoinfSNP_PAWR.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 75  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
}

locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])

# Load the genotype data
for (i in 1:length(chromosomes.to.analyze)) {
  chr <- chromosomes.to.analyze[i]  
  # Get chr data ---- 
  # read in position data for this chromosome
  position.file.name <- paste0(base.name, ".chr", chr, filename.text.middle, missing.genotype.threshold, filename.text.end, ".012.pos")
  pos.whole.chr <- read.table(position.file.name, col.names = c("chrom", "position"))
  # read in genotype data for this chromosome
  column_names <- c("null", paste("c", pos.whole.chr$chrom, pos.whole.chr$position, sep="."))
  genotype.file.name <- paste0(base.name, ".chr", chr, filename.text.middle, missing.genotype.threshold, filename.text.end, ".012NA")
  geno<-read.table(genotype.file.name, nrows = num.individuals, colClasses = "integer", col.names = column_names)
  loci_count <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
  # read in individual names for this chromosome dataset
  individuals.file.name <- paste0(base.name, ".chr", chr, filename.text.middle, missing.genotype.threshold, filename.text.end, ".012.indv")
  ind<-read.table(individuals.file.name)
}   


#set up an empty vector to hold the data
heterozygosities <- NULL
#For each row of genotypes, calculate heterozygosity (the number of "1" genotype calls divided by the total number of non-missing calls for that sample))
for(i in 1:nrow(geno)){
	heterozygosities$heterozygosity[i]<- sum(geno[i,]=="1", na.rm=T)/sum(!is.na(geno[i,]))
}
heterozygosities$sample <- ind
heterozygosities$location <- locations
heterozygosities <- as.data.frame(heterozygosities)

#compare heterozygosity of the different species
library(dplyr)
WIWR_het <- heterozygosities %>%
	filter(location.group=="WIWR")%>%
	summarise(mean(heterozygosity))
PAWR_het <- heterozygosities %>%
	filter(location.group=="PAWR")%>%
	summarise(mean(heterozygosity))


