# PCA Plotting

This file contains the code which was used to produce the PCA figures shown in the manuscript.

**Input**: `*.012NA`, `*.indv`, and `*.pos` files containing the genotype data for all samples ,produced in the first step. These files are provided in the folder `PAWR_WIWR_012NA_files` in this repository. This R code also uses functions that are provided in the R script `genomics_R_functions.R`, which is provided in this repository in the folder `2_PCA_data`. It also requires metadata files, which are provided in the folder `2_PCA_DATA` in this repository.

**Output:**: PCA figures and PCA data

Note that the R code, provided here, as well as other R code in this repository, will require a few changes before it can be run on a different machine: changing the path to the working directory, source function file, and input data files.

# PCA of whole genome data

First, here is the code used to produce a PCA of the whole genome data with all samples. It uses functions that are contained in the file `genomics_R_functions.R`, which is provided in this repository in the folder `2_PCA_data`.

This dataset uses only variable sites, with singletons removed (as these will not be useful for PCA of relationships between samples).

Note that this code has several lines which are not used to generate the PCA plots, but contain options used for other figures/analyses. 

```R
#here is the script I will use for the whole genome PCA

# Introduction ----
# PAWR_WIWR_GBS_R_analysis_script.R
# Started by Darren Irwin on 5 Feb 2018.
# Edited for whole SNP set PCA by Else Mikkelsen Oct 8 2018
# Based on R code written for Greenish Warbler analysis, and then the NA warbler analyses.

# Initial setup ----

# set directory where files stored
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")

# Load functions
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/genomics_R_functions.R")

##########genome wide PCA, PAWR vs WIWR############

#this is a PCA using a dataset of SNPs with singletons removed, with all chromosomes except for the sex chromosomes

# choose the chromosomes to analyze in this run
chromosomes.to.analyze <- c("autoinfSNP")  #autoinfSNP is a dataset containing all autsomes concatenated together, with only "informative" sites (no invariant sites)

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
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors <- c("blue", "red", "purple", "grey")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR")
  group.colors.WC84_Fst <- c("purple")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR")
  group.colors.pi <- c("blue", "red")
}

# Option to focus on a region of chromosome ----
focus.region <- F  # choose T for a subset of the chromosome, F for the whole thing)
if (focus.region==T){
  position.min <-  1500000
  position.max <- 1750000
}

calculate_or_load_stats <- 3  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) load per-site and windowed data from file

# Option to save the per-site stats
saveSiteInfo <- F # TRUE     # TRUE   # If TRUE, will save a file for per-site stats

saveWindowedStats <- F # TRUE     #TRUE   # If TRUE, will save a file for per-window stats

load.rolling.means <- T   #FALSE     #   FALSE    # If TRUE, will load rolling mean data (rather than calculate)

locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])


# MAIN LOOP ----
# -----
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
  
  # Get metadata ----
  
  ind_with_locations <- cbind(ind,locations) 
  print(ind_with_locations)    
  print("check first two columns to make sure the same")
  # combine metadata with genotype data
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  # If need to filter out individuals, based on low read number:
  filter <- T
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20),] #I am filtering this individual out because it is not independant of its close relative
  } else if (1==1) {
    combo.NApass.whole.chr <- combo
  }
  
  # Get region text ----
  if (focus.region==F) {
    position.min <- 1
    position.max <- pos.whole.chr$position[length(pos.whole.chr$position)]
    pos <- pos.whole.chr
    combo.NApass <- combo.NApass.whole.chr
    region.text <- paste0("Chr",chr,"_whole")
  } else if (focus.region==T) {
    selection <- pos.whole.chr$position >= position.min & pos.whole.chr$position <= position.max
    pos <- pos.whole.chr[selection,]
    selection <- c(rep(TRUE, times=num_loc_cols), selection)
    combo.NApass <- combo.NApass.whole.chr[, selection]
    region.text <- paste0("Chr",chr,"_from_",position.min,"_to_",position.max)
  }
  
  # Make site stats ---- 
  if (calculate_or_load_stats==1) {
    # Calculate allele freqs and sample sizes (use column Fst_group)
    temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
    freqs <- temp.list$freqs
    sample_size <- temp.list$sample_size
    rm(temp.list)
    print("Calculated population allele frequencies and sample sizes")
    
    # calculate nucleotide diversity (pi) at each site for each population
    site_pi <- getSitePi(freqs) 
    print("Calculated population pi values")
    
    # calculate Dxy at each site, between pairs of groups
    Dxy <- getDxy(freqs, groups)
    print("Calculated Dxy values")
    
    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    temp.list <- getWC84Fst(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
    WC84_Fst <- temp.list$WC84_Fst
    WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
    WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
    rm(temp.list)
    print("Calculated WC84_Fst values")
    
    if (saveSiteInfo == TRUE) {  # save the per-site stats, if chosen to in Intro section
      save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
      print("Saved summary site stats")
      WC84_Fst_among <- WC84_Fst[rownames(WC84_Fst)=="Fst_among",]
      print(paste0(length(WC84_Fst_among)," markers in total (",sum(!is.na(WC84_Fst_among))," variable and ",sum(is.na(WC84_Fst_among)), " invariant)"))
    } else print("Site stats not saved")
    
  } else if (calculate_or_load_stats==2 | calculate_or_load_stats==3) {
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print("Loaded saved summary stats")
  }
  # Make windowed stats ---- 
  if (calculate_or_load_stats==1 | calculate_or_load_stats==2) {
    
    # calculate windowed pi, in whole windows starting on left side of chromosome
    temp.list <- getWindowedPi(site_pi, pos, window_size, step_size)
    rolling.mean.pos.pi <- temp.list$rolling.mean.pos.pi
    rolling.mean.pi <- temp.list$rolling.mean.pi
    rm(temp.list)
    
    # calculate windowed Dxy, in whole windows starting on left side of chromosome
    temp.list <- getWindowedDxy(Dxy, pos, window_size, step_size)
    rolling.mean.pos.Dxy <- temp.list$rolling.mean.pos.Dxy
    rolling.mean.Dxy <- temp.list$rolling.mean.Dxy
    rm(temp.list)
    
    # calculate windowed Fst according to according to Weir&Cockerham1984 
    # (with sample size and pop number correction),
    # calculated as windowed numerator over windowed denominator.
    temp.list <- getWindowedWC84_Fst(WC84_Fst_numerator, WC84_Fst_denominator, pos, window_size, step_size)
    rolling.mean.pos.WC84_Fst <- temp.list$rolling.mean.pos.WC84_Fst
    rolling.mean.WC84_Fst <- temp.list$rolling.mean.WC84_Fst
    rm(temp.list)
  }
  
  # Save or load ----
  # save the rolling mean data, if chosen in Intro section
  if (saveWindowedStats == TRUE) {
    save(rolling.mean.pos.pi, rolling.mean.pi, rolling.mean.pos.Dxy, rolling.mean.Dxy, rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst, file=paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
    print("Saved rolling mean stats")
  }
  
  # load the rolling mean data, if chosen in Intro section:
  if (load.rolling.means == TRUE) {
    load(paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
  }
  
  # Make plot for chr ----
  # make ggplots for quick inspection of rolling mean results
  makeRollingMeanPlots(rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst,
                       rolling.mean.pos.Dxy, rolling.mean.Dxy,
                       rolling.mean.pos.pi, rolling.mean.pi, 
                       group.colors.pi, groups.to.plot.pi,
                       group.colors.Dxy, groups.to.plot.Dxy,
                       group.colors.WC84_Fst, groups.to.plot.WC84_Fst,
                       region.text)
  
}   
# End main loop ----


###This section performs a PCA plot for pacificus and hiemalis#####
## Make PCA plot ----
# choose only loci that are variable in the dataset (SNPs), and (optionally) above an Fst threshhold
# groups and colors determined in Intro section, under groups.to.plot.PCA and group.colors.PCA
Fst.filter <- F #F turns off the Fst filter
Fst.cutoff <- 0.01
# choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Fst_among"
#groups.to.compare <- "nigrifrons_auduboni"
axes <- 3

#select the colours to plot each population
groups.to.plot.PCA <- c("PAWR", "WIWR", "Hybrid")
group.colors.PCA <- c("blue", "red", "purple")

#store the PCA results into a variable
PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text, groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=F)

# PCA_results is a list containing var_explained, scores, and data 
PCA_results$var_explained

#to check PC3
#quartz()
#plot(PCA_results$scores[,1], PCA_results$scores[,3],  asp=1, pch=23, cex=2, xlab="PC1", ylab="PC3")
```

# PCA of T. pacificus samples

Next, here is the code for making the figure of PCA for *T. pacificus* without *T. hiemalis*. The goal was to determine if the samples from the range of *T. p. salebrosus* can be distinguished from the *T. p. pacificus* samples.

```R
#here is the script I will use for the whole genome PCA

# Introduction ----
# PAWR_WIWR_GBS_R_analysis_script.R
# Started by Darren Irwin on 5 Feb 2018.
# Edited for whole SNP set PCA by Else Mikkelsen Oct 8 2018
# Based on R code written for Greenish Warbler analysis, and then the NA warbler analyses.

# Initial setup ----

# set directory where files stored
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")

# Load functions
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/genomics_R_functions.R") #differs from old code at line 205

###Now this code makes a PCA of PAWR pop structure#####
#it will colour the points by subspecies

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
  metadata.file <- "PAWRWIWR_Metadata_autoinfSNP_PAWR.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 75  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("pacificus", "WIWR", "Hybrid", "salebrosus")
  group.colors <- c("blue", "red", "purple", "cadetblue1")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("pacificus_WIWR", "pacificus_Hybrid", "WIWR_Hybrid")
  group.colors.WC84_Fst <- c("purple", "blue", "red")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("pacificus", "WIWR", "Hybrid")
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
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors <- c("blue", "red", "purple", "grey")
  #groups <- c("PAWR", "WIWR", "Hybrid")
  #group.colors <- c("blue", "red", "grey")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  #groups.to.plot.WC84_Fst <- c("PAWR_WIWR", "PAWR_Hybrid", "WIWR_Hybrid")
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR", "PAWR_Hybrid", "PAWR_MAWR")
  group.colors.WC84_Fst <- c("purple", "blue", "grey")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR")
  group.colors.pi <- c("blue", "red")
}

# Option to focus on a region of chromosome ----
focus.region <- F  # choose T for a subset of the chromosome, F for the whole thing)
if (focus.region==T){
  position.min <-  1500000
  position.max <- 1750000
}

calculate_or_load_stats <- 3  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) load per-site and windowed data from file

# Option to save the per-site stats
saveSiteInfo <- F # TRUE     # TRUE   # If TRUE, will save a file for per-site stats

saveWindowedStats <- F # TRUE     #TRUE   # If TRUE, will save a file for per-window stats

load.rolling.means <- T   #FALSE     #   FALSE    # If TRUE, will load rolling mean data (rather than calculate)

locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])

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
  
  # Get metadata ----
  
  ind_with_locations <- cbind(ind,locations) 
  print(ind_with_locations)    
  print("check first two columns to make sure the same")
  # combine metadata with genotype data
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  # If need to filter out individuals, based on low read number:
  filter <- T
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20),] #I am filtering this individual out because it is not independant of its close relative
  } else if (1==1) {
    combo.NApass.whole.chr <- combo
  }
  
  # Get region text ----
  if (focus.region==F) {
    position.min <- 1
    position.max <- pos.whole.chr$position[length(pos.whole.chr$position)]
    pos <- pos.whole.chr
    combo.NApass <- combo.NApass.whole.chr
    region.text <- paste0("Chr",chr,"_whole")
  } else if (focus.region==T) {
    selection <- pos.whole.chr$position >= position.min & pos.whole.chr$position <= position.max
    pos <- pos.whole.chr[selection,]
    selection <- c(rep(TRUE, times=num_loc_cols), selection)
    combo.NApass <- combo.NApass.whole.chr[, selection]
    region.text <- paste0("Chr",chr,"_from_",position.min,"_to_",position.max)
  }
  
  # Make site stats ---- 
  if (calculate_or_load_stats==1) {
    # Calculate allele freqs and sample sizes (use column Fst_group)
    temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
    freqs <- temp.list$freqs
    sample_size <- temp.list$sample_size
    rm(temp.list)
    print("Calculated population allele frequencies and sample sizes")
    
    # calculate nucleotide diversity (pi) at each site for each population
    site_pi <- getSitePi(freqs) 
    print("Calculated population pi values")
    
    # calculate rownames for pairwise comparisons, for use in Dxy and Fst matrices:
    # rownames <- getPairwiseNames(groups)   # NOT NEEDED HERE since called from getDxy
    
    # calculate Dxy at each site, between pairs of groups
    Dxy <- getDxy(freqs, groups)
    print("Calculated Dxy values")
    
    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    temp.list <- getWC84Fst(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
    WC84_Fst <- temp.list$WC84_Fst
    WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
    WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
    rm(temp.list)
    print("Calculated WC84_Fst values")
    
    if (saveSiteInfo == TRUE) {  # save the per-site stats, if chosen to in Intro section
      save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
      print("Saved summary site stats")
      WC84_Fst_among <- WC84_Fst[rownames(WC84_Fst)=="Fst_among",]
      print(paste0(length(WC84_Fst_among)," markers in total (",sum(!is.na(WC84_Fst_among))," variable and ",sum(is.na(WC84_Fst_among)), " invariant)"))
    } else print("Site stats not saved")
    
  } else if (calculate_or_load_stats==2 | calculate_or_load_stats==3) {
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print("Loaded saved summary stats")
  }
  
  # Make windowed stats ---- 
  if (calculate_or_load_stats==1 | calculate_or_load_stats==2) {
    
    # calculate windowed pi, in whole windows starting on left side of chromosome
    temp.list <- getWindowedPi(site_pi, pos, window_size, step_size)
    rolling.mean.pos.pi <- temp.list$rolling.mean.pos.pi
    rolling.mean.pi <- temp.list$rolling.mean.pi
    rm(temp.list)
    
    # calculate windowed Dxy, in whole windows starting on left side of chromosome
    temp.list <- getWindowedDxy(Dxy, pos, window_size, step_size)
    rolling.mean.pos.Dxy <- temp.list$rolling.mean.pos.Dxy
    rolling.mean.Dxy <- temp.list$rolling.mean.Dxy
    rm(temp.list)
    
    # calculate windowed Fst according to according to Weir&Cockerham1984 
    # (with sample size and pop number correction),
    # calculated as windowed numerator over windowed denominator.
    temp.list <- getWindowedWC84_Fst(WC84_Fst_numerator, WC84_Fst_denominator, pos, window_size, step_size)
    rolling.mean.pos.WC84_Fst <- temp.list$rolling.mean.pos.WC84_Fst
    rolling.mean.WC84_Fst <- temp.list$rolling.mean.WC84_Fst
    rm(temp.list)
  }
  
  # Save or load ----
  # save the rolling mean data, if chosen in Intro section
  if (saveWindowedStats == TRUE) {
    save(rolling.mean.pos.pi, rolling.mean.pi, rolling.mean.pos.Dxy, rolling.mean.Dxy, rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst, file=paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
    print("Saved rolling mean stats")
  }
  
  # load the rolling mean data, if chosen in Intro section:
  if (load.rolling.means == TRUE) {
    load(paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
  }
  
  # Make plot for chr ----
  # make ggplots for quick inspection of rolling mean results
  makeRollingMeanPlots(rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst,
                       rolling.mean.pos.Dxy, rolling.mean.Dxy,
                       rolling.mean.pos.pi, rolling.mean.pi, 
                       group.colors.pi, groups.to.plot.pi,
                       group.colors.Dxy, groups.to.plot.Dxy,
                       group.colors.WC84_Fst, groups.to.plot.WC84_Fst,
                       region.text)
  
}   
# End main loop ----


## Make PCA plot ----
# choose only loci that are variable in the dataset (SNPs), and (optionally) above an Fst threshhold
# groups and colors determined in Intro section, under groups.to.plot.PCA and group.colors.PCA
Fst.filter <- F
Fst.cutoff <- 0.01
# choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Fst_among"
#groups.to.compare <- "nigrifrons_auduboni"
axes <- 3

groups.to.plot.PCA <- c("pacificus", "salebrosus")
group.colors.PCA <- c("blue", "cadetblue1")

PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                       groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=F)

# PCA_results is a list containing var_explained, scores, and data 
PCA_results$var_explained

PCA_results$scores[,1:3]

#to check PC3
#quartz()
#plot(PCA_results$scores[,1], PCA_results$scores[,3],  asp=1, pch=23, cex=2, xlab="PC1", ylab="PC3")
```

Before leaving the question of *T. pacificus salebrosus* and *T. pacificus pacificus*, I would like to see where the most highly differentiated SNP is.
```R
# Introduction ----
# PAWR_WIWR_GBS_R_analysis_script.R
# Started by Darren Irwin on 5 Feb 2018.
# Edited for whole SNP set PCA by Else Mikkelsen Oct 8 2018
# Based on R code written for Greenish Warbler analysis, and then the NA warbler analyses.

# Initial setup ----

# set directory where files stored
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")

# Load functions
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/genomics_R_functions.R") #differs from old code at line 205

###Now this code makes a PCA of PAWR pop structure#####
#it will colour the points by subspecies

chromosomes.to.analyze <- c("20")  

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
  metadata.file <- "PAWRWIWR_Metadata_subspecies.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 77  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 1
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("Pacificus", "WIWR", "Hybrid", "Salebrosus")
  group.colors <- c("blue", "red", "purple", "cadetblue1")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("Pacificus_Salebrosus")
  group.colors.WC84_Fst <- c("blue")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("Pacificus", "Salebrosus")
  group.colors.pi <- c("blue", "cadetblue1")
  
}

# Option to focus on a region of chromosome ----
focus.region <- T  # choose T for a subset of the chromosome, F for the whole thing)
if (focus.region==T){
  position.min <- 13700000
  position.max <- 13900000
}

calculate_or_load_stats <- 1  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) load per-site and windowed data from file

# Option to save the per-site stats
saveSiteInfo <- F # TRUE     # TRUE   # If TRUE, will save a file for per-site stats

saveWindowedStats <- F # TRUE     #TRUE   # If TRUE, will save a file for per-window stats

load.rolling.means <- F   #FALSE     #   FALSE    # If TRUE, will load rolling mean data (rather than calculate)


locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])

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
  
  # Get metadata ----
  
  ind_with_locations <- cbind(ind,locations) 
  print(ind_with_locations)    
  print("check first two columns to make sure the same")
  # combine metadata with genotype data
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  # If need to filter out individuals, based on low read number:
  filter <- T
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20),] #I am filtering this individual out because it is not independant of its close relative
  } else if (1==1) {
    combo.NApass.whole.chr <- combo
  }
  
  # Get region text ----
  if (focus.region==F) {
    position.min <- 1
    position.max <- pos.whole.chr$position[length(pos.whole.chr$position)]
    pos <- pos.whole.chr
    combo.NApass <- combo.NApass.whole.chr
    region.text <- paste0("Chr",chr,"_whole")
  } else if (focus.region==T) {
    selection <- pos.whole.chr$position >= position.min & pos.whole.chr$position <= position.max
    pos <- pos.whole.chr[selection,]
    selection <- c(rep(TRUE, times=num_loc_cols), selection)
    combo.NApass <- combo.NApass.whole.chr[, selection]
    region.text <- paste0("Chr",chr,"_from_",position.min,"_to_",position.max)
  }
  
  # Make site stats ---- 
  if (calculate_or_load_stats==1) {
    # Calculate allele freqs and sample sizes (use column Fst_group)
    temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
    freqs <- temp.list$freqs
    sample_size <- temp.list$sample_size
    rm(temp.list)
    print("Calculated population allele frequencies and sample sizes")
    
    # calculate nucleotide diversity (pi) at each site for each population
    site_pi <- getSitePi(freqs) 
    print("Calculated population pi values")
    
    # calculate rownames for pairwise comparisons, for use in Dxy and Fst matrices:
    # rownames <- getPairwiseNames(groups)   # NOT NEEDED HERE since called from getDxy
    
    # calculate Dxy at each site, between pairs of groups
    Dxy <- getDxy(freqs, groups)
    print("Calculated Dxy values")
    
    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    temp.list <- getWC84Fst(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
    WC84_Fst <- temp.list$WC84_Fst
    WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
    WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
    rm(temp.list)
    print("Calculated WC84_Fst values")
    
    if (saveSiteInfo == TRUE) {  # save the per-site stats, if chosen to in Intro section
      save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
      print("Saved summary site stats")
      WC84_Fst_among <- WC84_Fst[rownames(WC84_Fst)=="Fst_among",]
      print(paste0(length(WC84_Fst_among)," markers in total (",sum(!is.na(WC84_Fst_among))," variable and ",sum(is.na(WC84_Fst_among)), " invariant)"))
    } else print("Site stats not saved")
    
  } else if (calculate_or_load_stats==2 | calculate_or_load_stats==3) {
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print("Loaded saved summary stats")
  }
  
  # Make windowed stats ---- 
  if (calculate_or_load_stats==1 | calculate_or_load_stats==2) {
    
    # calculate windowed pi, in whole windows starting on left side of chromosome
    temp.list <- getWindowedPi(site_pi, pos, window_size, step_size)
    rolling.mean.pos.pi <- temp.list$rolling.mean.pos.pi
    rolling.mean.pi <- temp.list$rolling.mean.pi
    rm(temp.list)
    
    # calculate windowed Dxy, in whole windows starting on left side of chromosome
    temp.list <- getWindowedDxy(Dxy, pos, window_size, step_size)
    rolling.mean.pos.Dxy <- temp.list$rolling.mean.pos.Dxy
    rolling.mean.Dxy <- temp.list$rolling.mean.Dxy
    rm(temp.list)
    
    # calculate windowed Fst according to according to Weir&Cockerham1984 
    # (with sample size and pop number correction),
    # calculated as windowed numerator over windowed denominator.
    temp.list <- getWindowedWC84_Fst(WC84_Fst_numerator, WC84_Fst_denominator, pos, window_size, step_size)
    rolling.mean.pos.WC84_Fst <- temp.list$rolling.mean.pos.WC84_Fst
    rolling.mean.WC84_Fst <- temp.list$rolling.mean.WC84_Fst
    rm(temp.list)
  }
  
  # Save or load ----
  # save the rolling mean data, if chosen in Intro section
  if (saveWindowedStats == TRUE) {
    save(rolling.mean.pos.pi, rolling.mean.pi, rolling.mean.pos.Dxy, rolling.mean.Dxy, rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst, file=paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
    print("Saved rolling mean stats")
  }
  
  # load the rolling mean data, if chosen in Intro section:
  if (load.rolling.means == TRUE) {
    load(paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
  }
  
  # Make plot for chr ----
  # make ggplots for quick inspection of rolling mean results
  makeRollingMeanPlots(rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst,
                       rolling.mean.pos.Dxy, rolling.mean.Dxy,
                       rolling.mean.pos.pi, rolling.mean.pi, 
                       group.colors.pi, groups.to.plot.pi,
                       group.colors.Dxy, groups.to.plot.Dxy,
                       group.colors.WC84_Fst, groups.to.plot.WC84_Fst,
                       region.text)
  
}   
# End main loop ----


```


# PCA of *T. hiemalis* samples

The goal is to determine whether there is any genetic structure amongst the *T. hiemalis* samples: particularly, whether the eastern samples will cluster apart from the western samples. One complication is that three eastern samples are female, and no western samples are female. Including any sex-linked markers would be a problem, as it could introduce a false signal of differentiation between the eastern and western samples. Because of this, it is important to remove the W-linked region mapping to chromosome 8. To be conservative, this script excluded the entirety of chromosome 8. This makes it possible to distinguish population structure from sex differences.

To run PCA on the *hiemalis* individuals, it is necessary to have a datafile that has no columns of entirely missing data. This was not a problem with the *pacificus* analysis, as no columns were completely missing for all *pacificus* samples. Unfortunately, there was at least one locus that was missing in all *hiemalis* samples. A new dataset was made by removing all *pacificus* (and the two hybrids) from the dataset, and then removing all invariant loci, and making new 012NA files from this dataset. This dataset is provided with the other 012NA files.

```R
# Introduction ----
# PAWR_WIWR_GBS_R_analysis_script.R
# Started by Darren Irwin on 5 Feb 2018.
# Edited for whole SNP set PCA by Else Mikkelsen Oct 8 2018
# Based on R code written for Greenish Warbler analysis, and then the NA warbler analyses.

# Initial setup ----

# set directory where files stored
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")

# Load functions
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/genomics_R_functions.R") #differs from old code at line 205


#chromosomes.to.analyze <- c("autoinfSNPWIWR")  #This would run with chr8 not removed, and PC1 would capture the variance between males and females due to the chr8 translocation
chromosomes.to.analyze <- c("autoinfSNPno8")  #removing chr8 removes the confounding variation created by the fact that 3 eastern samples are female and no western samples are female

Analysis_set <- 1  # 1: all samples, only SNPs;    2: all samples, with invariant sites
if (Analysis_set == 1) {
  # choose path and filename for the 012NA files
  base.name <- "PAWR_WIWR_012NA_files/PAWR_WIWR.genotypes.SNPs_only"
  filename.text.middle <- ".max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 30 
  filename.text.end <-".MQ20.lowHet.tab"
  #tag.name <- ".no8_chrWIWR."   # choose a tag name for this analysis
  tag.name <- ".all_chr."   # choose a tag name for this analysis
  # indicate name of metadata file, a text file with these column headings:
  # ID  location  group Fst_group plot_order
  metadata.file <- "PAWRWIWR_Metadata_autoinfSNP_WIWR.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 14  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("WIWR", "WIWRE")
  group.colors <- c("red", "blue")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("WIWR_WIWRE")
  group.colors.WC84_Fst <- c("red")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("WIWR", "WIWRE")
  group.colors.pi <- c("red", "blue")
  
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
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors <- c("blue", "red", "purple", "grey")
  #groups <- c("PAWR", "WIWR", "Hybrid")
  #group.colors <- c("blue", "red", "grey")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  #groups.to.plot.WC84_Fst <- c("PAWR_WIWR", "PAWR_Hybrid", "WIWR_Hybrid")
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR", "PAWR_Hybrid", "PAWR_MAWR")
  group.colors.WC84_Fst <- c("purple", "blue", "grey")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR")
  group.colors.pi <- c("blue", "red")
  
}

# Option to focus on a region of chromosome ----
focus.region <- F  # choose T for a subset of the chromosome, F for the whole thing)
if (focus.region==T){
  position.min <-  1500000
  position.max <- 1750000
}

calculate_or_load_stats <- 1  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) load per-site and windowed data from file

# Option to save the per-site stats
saveSiteInfo <- F # TRUE     # TRUE   # If TRUE, will save a file for per-site stats

saveWindowedStats <- F # TRUE     #TRUE   # If TRUE, will save a file for per-window stats

load.rolling.means <- F   #FALSE     #   FALSE    # If TRUE, will load rolling mean data (rather than calculate)
  
locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])

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
  
  # Get metadata ----
  ind_with_locations <- cbind(ind,locations) 
  print(ind_with_locations)    
  print("check first two columns to make sure the same")
  # combine metadata with genotype data
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  # If need to filter out individuals, based on low read number:
  filter <- F
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20),] #I am filtering this individual out because it is not independant of its close relative
  } else if (1==1) {
    combo.NApass.whole.chr <- combo
  }
  
  # Get region text ----
  if (focus.region==F) {
    position.min <- 1
    position.max <- pos.whole.chr$position[length(pos.whole.chr$position)]
    pos <- pos.whole.chr
    combo.NApass <- combo.NApass.whole.chr
    region.text <- paste0("Chr",chr,"_whole")
  } else if (focus.region==T) {
    selection <- pos.whole.chr$position >= position.min & pos.whole.chr$position <= position.max
    pos <- pos.whole.chr[selection,]
    selection <- c(rep(TRUE, times=num_loc_cols), selection)
    combo.NApass <- combo.NApass.whole.chr[, selection]
    region.text <- paste0("Chr",chr,"_from_",position.min,"_to_",position.max)
  }
  
  # Make site stats ---- 
  if (calculate_or_load_stats==1) {
    # Calculate allele freqs and sample sizes (use column Fst_group)
    temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
    freqs <- temp.list$freqs
    sample_size <- temp.list$sample_size
    rm(temp.list)
    print("Calculated population allele frequencies and sample sizes")
    
    # calculate nucleotide diversity (pi) at each site for each population
    site_pi <- getSitePi(freqs) 
    print("Calculated population pi values")
    
    # calculate rownames for pairwise comparisons, for use in Dxy and Fst matrices:
    # rownames <- getPairwiseNames(groups)   # NOT NEEDED HERE since called from getDxy
    
    # calculate Dxy at each site, between pairs of groups
    Dxy <- getDxy(freqs, groups)
    print("Calculated Dxy values")
    
    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    temp.list <- getWC84Fst(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
    WC84_Fst <- temp.list$WC84_Fst
    WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
    WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
    rm(temp.list)
    print("Calculated WC84_Fst values")
    
      } 
}   
# End main loop ----

## Make PCA plot ----
# choose only loci that are variable in the dataset (SNPs), and (optionally) above an Fst threshhold
# groups and colors determined in Intro section, under groups.to.plot.PCA and group.colors.PCA
Fst.filter <- F
Fst.cutoff <- 0.01
# choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Fst_among"
#groups.to.compare <- "nigrifrons_auduboni"
axes <- 3

groups.to.plot.PCA <- c("WIWR", "WIWRE")
group.colors.PCA <- c("red", "blue")

PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                       groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=F)

# PCA_results is a list containing var_explained, scores, and data 
PCA_results$var_explained

PCA_results$scores[,1:3]

#to check PC3
#quartz()
#plot(PCA_results$scores[,1], PCA_results$scores[,3],  asp=1, pch=23, cex=2, xlab="PC1", ylab="PC3")
#plot(PCA_results$scores[,3], PCA_results$scores[,1], pch=23, cex=2, xlab="PC3", ylab="PC1", col=)
```

# Chromsome-by-chromosome PCA

Finally, this set of code will produce a PCA for each chromosome using all samples. This will allow investigation of each chromosome separately and could draw attention to any unique patterns. For example, this is how the W chromosome-linked duplication of chromosome 8 sequence was spotted.

```R
#here is the script I will use for the whole genome PCA

# Introduction ----
# PAWR_WIWR_GBS_R_analysis_script.R
# Started by Darren Irwin on 5 Feb 2018.
# Edited for whole SNP set PCA by Else Mikkelsen Oct 8 2018
# Based on R code written for Greenish Warbler analysis, and then the NA warbler analyses.

# Initial setup ----

# set directory where files stored
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")

# Load functions
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/genomics_R_functions.R") #differs from old code at line 205


#select which chromosomes to analyze
chromosomes.to.analyze <- c("Z") #to analyze the sex chromosome 
chromosomes.to.analyze <- c("1", "1A","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28") #to analyze the autosomes

Analysis_set <- 1  # 1: all samples, only SNPs;    2: all samples, with invariant sites

if (Analysis_set == 1) {
  # choose path and filename for the 012NA files
  base.name <- "PAWR_WIWR_012NA_files/PAWR_WIWR.genotypes.SNPs_only"
  filename.text.middle <- ".max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 30 
  filename.text.end <-".MQ20.lowHet.tab"
  tag.name <- ".all_samples."   # choose a tag name for this analysis
  # indicate name of metadata file, a text file with these column headings:
  # ID  location  group Fst_group plot_order
  metadata.file <- "PAWRWIWR_Metadata.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 77  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors <- c("blue", "red", "purple", "grey")
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
  metadata.file <- "PAWRWIWR_Metadata.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 75  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors <- c("blue", "red", "purple", "grey")
  group_count <- length(groups)
  
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR")
  group.colors.WC84_Fst <- c("purple", "blue", "grey")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR")
  group.colors.pi <- c("blue", "red")
}

# Option to focus on a region of chromosome ----
focus.region <- F  # choose T for a subset of the chromosome, F for the whole thing)

calculate_or_load_stats <- 3  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) load per-site and windowed data from file

# Option to save the per-site stats
saveSiteInfo <- F # TRUE     # TRUE   # If TRUE, will save a file for per-site stats

saveWindowedStats <- F # TRUE     #TRUE   # If TRUE, will save a file for per-window stats

load.rolling.means <- T   #FALSE     #   FALSE    # If TRUE, will load rolling mean data (rather than calculate)

locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])

# MAIN LOOP ----
# -----
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
  
  # Get metadata ----
  
  ind_with_locations <- cbind(ind,locations) 
  print(ind_with_locations)    
  print("check first two columns to make sure the same")
  # combine metadata with genotype data
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  # If need to filter out individuals, based on low read number:
  filter <- T
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20),] #I am filtering this individual out because it is not independant of its close relative
  } else if (1==1) {
    combo.NApass.whole.chr <- combo
  }
  
  # Get region text ----
  if (focus.region==F) {
    position.min <- 1
    position.max <- pos.whole.chr$position[length(pos.whole.chr$position)]
    pos <- pos.whole.chr
    combo.NApass <- combo.NApass.whole.chr
    region.text <- paste0("Chr",chr,"_whole")
  } else if (focus.region==T) {
    selection <- pos.whole.chr$position >= position.min & pos.whole.chr$position <= position.max
    pos <- pos.whole.chr[selection,]
    selection <- c(rep(TRUE, times=num_loc_cols), selection)
    combo.NApass <- combo.NApass.whole.chr[, selection]
    region.text <- paste0("Chr",chr,"_from_",position.min,"_to_",position.max)
  }
  
  # Make site stats ---- 
  if (calculate_or_load_stats==1) {
    # Calculate allele freqs and sample sizes (use column Fst_group)
    temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
    freqs <- temp.list$freqs
    sample_size <- temp.list$sample_size
    rm(temp.list)
    print("Calculated population allele frequencies and sample sizes")
    
    # calculate nucleotide diversity (pi) at each site for each population
    site_pi <- getSitePi(freqs) 
    print("Calculated population pi values")

    # calculate Dxy at each site, between pairs of groups
    Dxy <- getDxy(freqs, groups)
    print("Calculated Dxy values")
    
    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    temp.list <- getWC84Fst(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
    WC84_Fst <- temp.list$WC84_Fst
    WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
    WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
    rm(temp.list)
    print("Calculated WC84_Fst values")
    
    if (saveSiteInfo == TRUE) {  # save the per-site stats, if chosen to in Intro section
      save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
      print("Saved summary site stats")
      WC84_Fst_among <- WC84_Fst[rownames(WC84_Fst)=="Fst_among",]
      print(paste0(length(WC84_Fst_among)," markers in total (",sum(!is.na(WC84_Fst_among))," variable and ",sum(is.na(WC84_Fst_among)), " invariant)"))
    } else print("Site stats not saved")
    
  } else if (calculate_or_load_stats==2 | calculate_or_load_stats==3) {
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print("Loaded saved summary stats")
  }
  
  # Make windowed stats ---- 
  if (calculate_or_load_stats==1 | calculate_or_load_stats==2) {
    
    # calculate windowed pi, in whole windows starting on left side of chromosome
    temp.list <- getWindowedPi(site_pi, pos, window_size, step_size)
    rolling.mean.pos.pi <- temp.list$rolling.mean.pos.pi
    rolling.mean.pi <- temp.list$rolling.mean.pi
    rm(temp.list)
    
    # calculate windowed Dxy, in whole windows starting on left side of chromosome
    temp.list <- getWindowedDxy(Dxy, pos, window_size, step_size)
    rolling.mean.pos.Dxy <- temp.list$rolling.mean.pos.Dxy
    rolling.mean.Dxy <- temp.list$rolling.mean.Dxy
    rm(temp.list)
    
    # calculate windowed Fst according to according to Weir&Cockerham1984 
    # (with sample size and pop number correction),
    # calculated as windowed numerator over windowed denominator.
    temp.list <- getWindowedWC84_Fst(WC84_Fst_numerator, WC84_Fst_denominator, pos, window_size, step_size)
    rolling.mean.pos.WC84_Fst <- temp.list$rolling.mean.pos.WC84_Fst
    rolling.mean.WC84_Fst <- temp.list$rolling.mean.WC84_Fst
    rm(temp.list)
  }
  
  # Save or load ----
  # save the rolling mean data, if chosen in Intro section
  if (saveWindowedStats == TRUE) {
    save(rolling.mean.pos.pi, rolling.mean.pi, rolling.mean.pos.Dxy, rolling.mean.Dxy, rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst, file=paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
    print("Saved rolling mean stats")
  }
  
  # load the rolling mean data, if chosen in Intro section:
  if (load.rolling.means == TRUE) {
    load(paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
  }
  
  # Make PCA plot for chr ----
  groups.to.plot.PCA <- c("PAWR", "WIWR", "Hybrid")
  group.colors.PCA <- c("blue", "red", "purple")
  plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
          groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=F)
}   
# End main loop ----
```

Note on missing data:  
When imputing values for PCA, it is best that the samples do not contain too large a fraction of missing data.  
Here is how I checked the amount of missing data:  

```bash
cat PAWR_WIWR.genotypes.SNPs_only.chrautoinfSNP.max2allele_noindel.maxmiss30.MQ20.lowHet.tab.012NA | while read i ; do echo "$i" | grep -o "NA" | wc -l ; done
```
Output:
```
8101
13636
12019
15846
8843
9663
7648
10706
9834
9572
8937
9965
7773
9213
8340
9828
9241
7804
8004
12675
7045
6388
7380
8323
10180
9241
9014
9082
9622
7786
8186
7728
8109
8134
7988
8344
7831
8853
8001
7411
6817
7754
7238
7930
7017
7629
9584
9446
7261
6884
7980
7580
9215
9356
9045
7499
7840
7393
8302
8180
7677
7635
7371
8455
10945
7534
8874
8497
8080
7677
8383
8885
7156
8123
10874
```
All but two samples have less than 10% (of total 127,779 sites) missing data, and so SVDimpute should perform well.