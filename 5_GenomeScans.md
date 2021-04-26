# Genome scans

This file details the pipeline used to perform genome scans of Fst, pi_within, and pi_between for the Pacific and Winter Wren dataset.

**Input**: `*.012NA`, `*.indv`, and `*.pos` files containing the genotype data for all samples ,produced in the first step. These files are provided in the folder `2_PCA_DATA` in this repository. It also requires metadata files which are provided in the same folder. This R code also uses functions that are provided in the R script `genomics_R_functions.R`.

**Output**: estimates of Fst, pi_within, and pi_between, for windows across the genome.

# Main loop

The first step of this script runs through the data to calculate Fst, pi, and Dxy (pi_between). This takes a very long time to run, and so alternately, the data can be pre-loaded using the WindowStats_from_R.R files (the path to these files will need to be modified for the script to run on a different machine). The path to the working directory will need to be modified, as well as the path to the `Genomics_R_Functions.R` script

This script is set to pre-load the data using the line `calculate_or_load_stats <- 3`, this can be changed to re-calculate the stats from the data.

Once complete, this portion will create quick rolling mean plots for each chromosome seperately for inspection. I have placed each of these plots into the folder `5_GenomeScans_Plots` for reference.

```R
# Introduction ----
# PAWR_WIWR_GBS_R_analysis_script.R
# Started by Darren Irwin on 5 Feb 2018.
# Based on R code written for Greenish Warbler analysis, and then the NA warbler analyses.

# Initial setup ----

# set directory where files stored
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")

# Load functions
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/Genomics_R_Functions_changeddxytopib.R")

# choose the chromosomes to analyze in this run
#chromosomes.to.analyze <- c("24", "25", "26")  ### I have put just an example of a set of chromosomes here--I have already run these through the main loop below. 
#### To do them all at once (which would take a long time), use this line:
chromosomes.to.analyze <- c("1", "1A", "2", "3","4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15","17","18","19","20","21","22","23","24","25","26","27","28","Z")

Analysis_set <- 2  # 1: all samples, only SNPs;    2: all samples, with invariant sites

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
  # ID	location	group	Fst_group	plot_order
  metadata.file <- "PAWRWIWR_Metadata.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 77  # specify number of individuals in file (good for error checking)
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
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR", "PAWR_Hybrid", "WIWR_Hybrid")
  group.colors.WC84_Fst <- c("purple", "blue", "red")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR", "Hybrid")
  group.colors.pi <- c("blue", "red", "grey")
  
} else if (Analysis_set == 2) {
  # choose path and filename for the 012NA files
  base.name <- "PAWR_WIWR_012NA_files/PAWR_WIWR.genotypes.allSites"
  filename.text.middle <- ".infoSites.max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 30 
  filename.text.end <-".MQ20.lowHet.tab"
  tag.name <- ".all_samples_with_invariant_fixingpichangedxy."   # choose a tag name for this analysis
  # indicate name of metadata file, a text file with these column headings:
  # ID	location	group	Fst_group	plot_order
  metadata.file <- "PAWRWIWR_Metadata.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 77  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 10000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors <- c("blue", "red", "purple", "grey")
  #groups <- c("PAWR", "WIWR", "Hybrid")
  #group.colors <- c("blue", "red", "grey")
  group_count <- length(groups)
  
  
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR", "PAWR_Hybrid", "WIWR_Hybrid", "PAWR_MAWR", "WIWR_MAWR")
  #groups.to.plot.WC84_Fst <- c("PAWR_WIWR", "PAWR_Hybrid", "PAWR_MAWR")
  group.colors.WC84_Fst <- c("purple", "blue", "red", "grey", "orange")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors.pi <- c("blue", "red", "purple", "grey")
  
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
  filter <- F
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20,-163),]
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
    #site_pi <- getSitePi(freqs) #for uncorrected pi, biased low at low sample size
    site_pi_nb <- getSitePi_nb(freqs, sample_size) #corrected pi
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
      #save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R")) #uncorrected pi
      save(pos, freqs, sample_size, site_pi_nb, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
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
    #temp.list <- getWindowedPi(site_pi, pos, window_size, step_size) #uncorrected pi
    temp.list <- getWindowedPi(site_pi_nb, pos, window_size, step_size)
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

After reading in the data, we can produce a genome scan of Fst, pi_within, and pi_between combining all chromosomes in the genome. This will produce a plot that I have placed into the `5_GenomeScans_Plots` folder: `WC84_Fst_pib_piw_allchromosomes.pdf`
```R
# GENOME-WIDE plots -------
# read files for each chromosome, and plot Fst, Dxy, and pi for each chromosome in a single window

# specify groups for plotting, and their colors
groups.to.plot.WC84_Fst <- c("PAWR_WIWR") #We only have two main populations, so only one Fst scan to plot
group.colors.WC84_Fst <- c("purple")

groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
groups.to.plot.pi <- c("PAWR", "WIWR") #we have two main populations to plot seperately
group.colors.pi <- c("blue", "red")

#chromosomes.to.plot <- c("24", "25", "26")   # You can add more if you have processed more above.
chromosomes.to.plot <- c("1","1A","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","Z") #plots all the chromosomes
max_Dxy_axis <- 0.01  # height of Dxy axis
max_pi_axis <- 0.01  # height of pi axis
#transparency <- 0.2  # transparency of the colored area under the lines in 2016 code
transparency_Fst_Dxy <- 0.4 #2018 code
transparency_pi <- 0.2 #2018 code
line_transparency <- 0.8  # transparency of the colored lines
# plot Genome-wide plot of windowed Fst, Dxy, pi
plotGenomeFstDxyPi(base.name, tag.name, window_size, chromosomes.to.plot, 
                   max_Dxy_axis, max_pi_axis, transparency_Fst_Dxy, transparency_pi, line_transparency,
                   groups.to.plot.WC84_Fst, group.colors.WC84_Fst,
                   groups.to.plot.Dxy, group.colors.Dxy,
                   groups.to.plot.pi, group.colors.pi) #Note the newer 2018 version of the script has 2 transparency variables instead of 1, if using 2016 functions only use 1 transparency variable
``` 

Next, we can produce a figure comparing Fst, pi_between, and pi_within at each genomic window using a scatterplot. We can also test the correlation between pi_w and pi_between First, we compile data from all the autosomes, and then create the plot using this data.
```R
# Compile genome-wide info  ----
# compile windowed WC84_Fst, Dxy, pi for a bunch of chromosomes
# compile autosomal-genome-wide info:
chromosomes.to.combine <- c("1","1A","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28")
autosome.genome.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window.size)


# Dxy_MeanPi plot ----
# Scatterplot of Dxy vs. mean_pi using compiled genome-wide info .
# Color the points based on Fst.
group1 <- "PAWR"
group2 <- "WIWR"
max.x <- 0.010   
max.y <- 0.010   
cor.method <- "spearman"

group1 <- "PAWR"
group2 <- "WIWR"
test_correlation <- plotDxy_MeanPi(group1, group2, cor.method,
                       autosome.genome.rolling.stats,
                       max.x, max.y, FALSE, c("grey", "blue"))
test_correlation   # print test of correlation between pi_w and pi_b

#results:
#	Spearman's rank correlation rho
#
#data:  Dxy and mean_pi
#S = 36040454, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#0.86182 

plotFst_MeanPi(group1, group2, cor.method,
                       autosome.genome.rolling.stats,
                       0.01, 1, FALSE, c("grey", "blue"))
plotDxy_Fst(group1, group2, cor.method,
                       autosome.genome.rolling.stats,
                       0.01, 1, FALSE, c("grey", "blue"))

plotDxy_MeanPi<- function(group1, group2, cor.method,
                           autosome.genome.rolling.stats, 
                           max.x, max.y, color_bar, color_extremes) {
  groups.for.graph <- paste0(group1, "_", group2)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  mean_pi <- (pi_1 + pi_2) / 2
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
  Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  WC84_Fst.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == groups.for.graph)
  WC84_Fst.vector <- autosome.genome.rolling.stats$WC84_Fst[WC84_Fst.row.choice,]
  # color points as a gradient according to WC84_Fst:
  color.palette <- colorRampPalette(color_extremes)(100)  # makes a color scale with 100 colors
  max.WC84_Fst <- max(is.finite(WC84_Fst.vector))  # the next line will cause the max.WC84_Fst to be set to the top color
  color_scale <- (round(WC84_Fst.vector * 100) / max.WC84_Fst)
  color_scale[color_scale<1] <- 1
  colors <- color.palette[color_scale]  # makes a list of colors according to Fst values of the windows
  quartz(title=paste0("Scatterplot of windowed mean pi vs. Dxy between ", group1,"_", group2, sep=""), width=5, height=5.5)
  plot(x=Dxy, y=mean_pi, col=colors, xlim=c(0,max(c(max(Dxy)*1.04,max.x))), ylim=c(0, max(c(max(mean_pi)*1.04,max.y))),
       xaxp = c(0, max.x, 2), yaxp=c(0, max.y, 2), pch=19, cex=0.5, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  lines(c(0,1), c(0,1))
  label.size <- 1.4
  if (color_bar == T) {
    colorbar.plot(x=0, y=0.92*max.y, seq(1,100,by=1), col = color.palette, strip.width=0.025, strip.length=0.4, horizontal=T, adj.x = 0, adj.y = 0)
    text(x=0.22*max.x, y=0.93*max.y, bquote(italic('F')[ST]* " scale"), pos=3, cex=0.8)
    text(x=0.005*max.x, y=0.925*max.y, "0", pos=1, cex=0.7)
    text(x=0.415*max.x, y=0.925*max.y, "1", pos=1, cex=0.7)
  }
  title(xlab=expression("Between-group differentiation (" * italic(pi)[B] * ")"), line=2.5, cex.lab=label.size)
  title(ylab=expression("Mean within-group variation (" * italic(pi)[W] * ")"), line=2, cex.lab=label.size)
  test <- cor.test(Dxy, mean_pi, method=cor.method)
  return(test)
}
plotFst_MeanPi<- function(group1, group2, cor.method,
                           autosome.genome.rolling.stats, 
                           max.x, max.y, color_bar, color_extremes) {
  groups.for.graph <- paste0(group1, "_", group2)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  mean_pi <- (pi_1 + pi_2) / 2
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
  Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  WC84_Fst.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == groups.for.graph)
  WC84_Fst.vector <- autosome.genome.rolling.stats$WC84_Fst[WC84_Fst.row.choice,]
  # color points as a gradient according to WC84_Fst:
  color.palette <- colorRampPalette(color_extremes)(100)  # makes a color scale with 100 colors
  max.Dxy <- max(is.finite(Dxy))  # the next line will cause the max.WC84_Fst to be set to the top color
  color_scale <- round(Dxy * 10000)
  #color_scale[color_scale<1] <- 1
  colors <- color.palette[color_scale]  # makes a list of colors according to Fst values of the windows
  quartz(title=paste0("Scatterplot of windowed mean pi vs. Dxy between ", group1,"_", group2, sep=""), width=5, height=5.5)
  plot(x=mean_pi, y=WC84_Fst.vector, col=colors, xlim=c(0,max(c(max(Dxy)*1.04,max.x))), ylim=c(0, max(c(max(mean_pi)*1.04,max.y))),
       xaxp = c(0, max.x, 2), yaxp=c(0, max.y, 2), pch=19, cex=0.5, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  label.size <- 1.4
  if (color_bar == T) {
    colorbar.plot(x=0, y=0.92*max.y, seq(1,100,by=1), col = color.palette, strip.width=0.025, strip.length=0.4, horizontal=T, adj.x = 0, adj.y = 0)
    text(x=0.22*max.x, y=0.93*max.y, bquote(italic('F')[ST]* " scale"), pos=3, cex=0.8)
    text(x=0.005*max.x, y=0.925*max.y, "0", pos=1, cex=0.7)
    text(x=0.415*max.x, y=0.925*max.y, "1", pos=1, cex=0.7)
  }
  title(xlab=expression("Mean within-group variation (" * italic(pi)[W] * ")"), line=2, cex.lab=label.size)
  title(ylab=expression("Relative Differentiation (" * italic(F)[ST] * ")"), line=2, cex.lab=label.size)
}
plotDxy_Fst<- function(group1, group2, cor.method,
                           autosome.genome.rolling.stats, 
                           max.x, max.y, color_bar, color_extremes) {
  groups.for.graph <- paste0(group1, "_", group2)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  mean_pi <- (pi_1 + pi_2) / 2
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
  Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  WC84_Fst.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == groups.for.graph)
  WC84_Fst.vector <- autosome.genome.rolling.stats$WC84_Fst[WC84_Fst.row.choice,]
  # color points as a gradient according to WC84_Fst:
  color.palette <- colorRampPalette(color_extremes)(100)  # makes a color scale with 100 colors
  color_scale <- round(mean_pi * 10000) # / max.WC84_Fst)
  #color_scale[color_scale<1] <- 1
  colors <- color.palette[color_scale]  # makes a list of colors according to Fst values of the windows
  quartz(title=paste0("Scatterplot of windowed mean pi vs. Dxy between ", group1,"_", group2, sep=""), width=5, height=5.5)
  plot(x=Dxy, y=WC84_Fst.vector, col=colors, xlim=c(0,max(c(max(Dxy)*1.04,max.x))), ylim=c(0, max(c(max(mean_pi)*1.04,max.y))),
       xaxp = c(0, max.x, 2), yaxp=c(0, max.y, 2), pch=19, cex=0.5, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  label.size <- 1.4
  if (color_bar == T) {
    colorbar.plot(x=0, y=0.92*max.y, seq(1,100,by=1), col = color.palette, strip.width=0.025, strip.length=0.4, horizontal=T, adj.x = 0, adj.y = 0)
    text(x=0.22*max.x, y=0.93*max.y, bquote(italic('F')[ST]* " scale"), pos=3, cex=0.8)
    text(x=0.005*max.x, y=0.925*max.y, "0", pos=1, cex=0.7)
    text(x=0.415*max.x, y=0.925*max.y, "1", pos=1, cex=0.7)
  }
  title(xlab=expression("Between-group differentiation (" * italic(pi)[B] * ")"), line=2.5, cex.lab=label.size)
  title(ylab=expression("Relative Differentiation (" * italic(F)[ST] * ")"), line=2, cex.lab=label.size)
}
```

Next, we can run a boundary analysis to determine whether there is a difference in pi_within/pi_between above and below a given Fst threshold. I tested different boundaries in order to identify that Fst value for which 95% of the windows were below the threshold (this turned out to be Fst=0.54) and the Fst value for which 99% f the windows were below the threshold (this turned out to be Fst=0.74).

```bash
# Fst boundary analysis ----
# make a histogram of Fst for one comparison
# and compare pi and Dxy above and below an Fst cutoff
group1 <- "PAWR"         
group2 <- "WIWR" 
Fst_boundary <- 0 #this could be to examine across the whole genome
#Fst_boundary <- 0.54 #this could be the cutoff for the 5% most differentiated windows (eg, 0.91 for PAWR vs MAWR, 0.54 for WIWR vs PAWR)
#Fst_boundary <- 0.74 #this could be to examine at the cutoff for the 1% most differentiated windows 
Fst_boundary_stats <- Fst_hist_and_boundary_stats(group1, group2, 
                                                  autosome.genome.rolling.stats, 
                                                  Fst_boundary)
Fst_boundary_stats  # prints some useful statistics

#results for boundary= 0.54:
#$fraction_windows_above_Fst_boundary
#[1] 0.04995693
#
#$mean_Dxy_above_Fst_boundary
#[1] 0.003654259
#
#$mean_Dxy_below_Fst_boundary
#[1] 0.003932947
#
#$mean_pi_above_Fst_boundary
#[1] 0.001338474
#
#$mean_pi_below_Fst_boundary
#[1] 0.003072358


#results for boundary= 0.74:
#$fraction_windows_above_Fst_boundary
#[1] 0.01033592
#
#$mean_Dxy_above_Fst_boundary
#[1] 0.00369615
#
#$mean_Dxy_below_Fst_boundary
#[1] 0.003921352
#
#$mean_pi_above_Fst_boundary
#[1] 0.0007852529
#
#$mean_pi_below_Fst_boundary
#[1] 0.00300872
#
#$mean_pi_1_above_Fst_boundary
#[1] 0.0006227417
#
#$mean_pi_1_below_Fst_boundary
#[1] 0.002661776
#
#$mean_pi_2_above_Fst_boundary
#[1] 0.000947764
#
#$mean_pi_2_below_Fst_boundary
#[1] 0.003355663


#results for boundary= 0:
#$fraction_windows_above_Fst_boundary
#[1] 1
#
#$mean_Dxy_above_Fst_boundary
#[1] 0.003919024
#
#$mean_Dxy_below_Fst_boundary
#[1] NaN
#
#$mean_pi_above_Fst_boundary
#[1] 0.002985738
#
#$mean_pi_below_Fst_boundary
#[1] NaN
#
#$mean_pi_1_above_Fst_boundary
#[1] 0.002640701
#
#$mean_pi_1_below_Fst_boundary
#[1] NaN
#
#$mean_pi_2_above_Fst_boundary
#[1] 0.003330776
#
#$mean_pi_2_below_Fst_boundary
#[1] NaN


#Here is the code for the function:
Fst_hist_and_boundary_stats <- function(group1, group2, 
                                        autosome.genome.rolling.stats, 
                                        Fst_boundary) {
  group.pair <- paste0(group1, "_", group2)
  quartz(title=paste0("Fst histogram, based on all autosomes, ", group.pair, sep=""), width=8, height=3)
  Fst.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == group.pair)
  Fst.vector <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice,]
  hist(Fst.vector)
  response <- NULL
  response$fraction_windows_above_Fst_boundary  <- sum(Fst.vector >= Fst_boundary) / length(Fst.vector)
  
  quartz(title=paste0("Dxy histogram, based on all autosomes, ", group.pair, sep=""), width=8, height=3)
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == group.pair)
  Dxy.vector <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice,]
  hist(Dxy.vector)
  response$mean_Dxy_above_Fst_boundary <- mean(Dxy.vector[Fst.vector >= Fst_boundary])
  response$mean_Dxy_below_Fst_boundary <- mean(Dxy.vector[Fst.vector < Fst_boundary])
  pi.row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  pi.row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  pi.vector1 <- as.vector(autosome.genome.rolling.stats$pi[pi.row.choice.1,])
  pi.vector2 <- as.vector(autosome.genome.rolling.stats$pi[pi.row.choice.2,])
  mean_pi <- (pi.vector1 + pi.vector2) / 2
  response$mean_pi_above_Fst_boundary <- mean(mean_pi[Fst.vector >= Fst_boundary])
  response$mean_pi_below_Fst_boundary <- mean(mean_pi[Fst.vector < Fst_boundary])
  response$mean_pi_1_above_Fst_boundary <- mean(pi.vector1[Fst.vector >= Fst_boundary])
  response$mean_pi_1_below_Fst_boundary <- mean(pi.vector1[Fst.vector < Fst_boundary])
  response$mean_pi_2_above_Fst_boundary <- mean(pi.vector2[Fst.vector >= Fst_boundary])
  response$mean_pi_2_below_Fst_boundary <- mean(pi.vector2[Fst.vector < Fst_boundary])
  
  quartz(title=paste0("mean_pi histogram grouped by Fst, based on all autosomes, ", group.pair, sep=""), width=8, height=5)
  bin_size <- 0.00025
  tick_label_gap <- 4
  xrange = seq(0, max(mean_pi)+bin_size, bin_size)
  hist1 = hist(mean_pi[Fst.vector >= Fst_boundary],breaks=xrange,plot=F)$counts
  hist2 = hist(mean_pi[Fst.vector < Fst_boundary],breaks=xrange,plot=F)$counts
  barplot(rbind(hist1,hist2),col=2:4,axes=FALSE,space=0,las=1)
  at_tick <- seq(0,length(hist1), by=tick_label_gap)
  tick_label <- at_tick * bin_size
  axis(side = 2, pos = -0.2, cex.axis=1)  # the y axis
  axis(side = 1, pos = -0.2, at = at_tick, labels = tick_label, cex.axis=1)
  
  return(response)
}
```

Next, we can calculate the mean and median values for pi, dxy, and fst for each species/comparison. I will calculate both with and without the Z chromosome, as I am curious to see how it affects the estimate.
```R
#I want to know average pi_w, pi_d, and Fst

#PAWR
mean(rolling.mean.pi[1,]) # Z, 0.00213
mean(autosome.genome.rolling.stats$pi[1,]) #no Z 0.00264
median(rolling.mean.pi[1,]) # Z 0.00212
median(autosome.genome.rolling.stats$pi[1,]) #no Z 0.00261
#including chrZ lowers the average of pi_w

#WIWR
mean(rolling.mean.pi[2,]) #0.00303
mean(autosome.genome.rolling.stats$pi[2,]) #no Z 0.00333
median(rolling.mean.pi[2,]) #0.00292
median(autosome.genome.rolling.stats$pi[2,]) #no Z 0.00333

#MAWR
mean(rolling.mean.pi[4,]) #0.00129
mean(autosome.genome.rolling.stats$pi[4,]) #no Z 0.00105
median(rolling.mean.pi[4,]) #0.000541
median(autosome.genome.rolling.stats$pi[4,]) #no Z 0.000690

#PAWR vs WIWR
median(rolling.mean.Dxy[1,]) #PAWR WIWR  Z 0.003638421
median(autosome.genome.rolling.stats$Dxy[1,]) #PAWR WIWR no Z 0.003822752
mean(rolling.mean.Dxy[1,]) #PAWR WIWR  Z 0.003958162
mean(autosome.genome.rolling.stats$Dxy[1,]) #PAWRWIWR no Z 0.003919024

median(rolling.mean.WC84_Fst[1,]) #PAWR WIWR  Z 0.3497863
median(autosome.genome.rolling.stats$WC84_Fst[1,]) #PAWR WIWR no Z 0.2179398
mean(rolling.mean.WC84_Fst[1,]) #PAWR WIWR  Z 0.3593871
mean(autosome.genome.rolling.stats$WC84_Fst[1,]) #PAWR WIWR no Z 0.251662

# PAWR vs MAWR
median(rolling.mean.Dxy[3,]) #PAWR MAWR 0.01253483
mean(rolling.mean.Dxy[3,]) #PAWR MAWR 0.01267135
median(autosome.genome.rolling.stats$Dxy[3,]) #PAWR MAWR no Z 0.01189639
mean(autosome.genome.rolling.stats$Dxy[3,]) #PAWR MAWR no Z 0.01213406

median(rolling.mean.WC84_Fst[3,]) #PAWR MAWR 0.8457175
mean(rolling.mean.WC84_Fst[3,]) #PAWR MAWR 0.8415678
median(autosome.genome.rolling.stats$WC84_Fst[3,]) #PAWR MAWR no Z 0.7981926
mean(autosome.genome.rolling.stats$WC84_Fst[3,]) #PAWR MAWR no Z 0.7975066

#WIWR vs MAWR
median(rolling.mean.Dxy[5,]) #WIWR MAWR 0.01252245
mean(rolling.mean.Dxy[5,]) #WIWR MAWR 0.01284918
median(autosome.genome.rolling.stats$Dxy[5,]) #MAWR WIWR no Z 0.01196223
mean(autosome.genome.rolling.stats$Dxy[5,]) #MAWR WIWR no Z 0.01218476

median(rolling.mean.WC84_Fst[5,]) #WIWR MAWR 0.7962195
mean(rolling.mean.WC84_Fst[5,]) #WIWR MAWR 0.7842912
median(autosome.genome.rolling.stats$WC84_Fst[5,]) #MAWR WIWR no Z 0.7544667
mean(autosome.genome.rolling.stats$WC84_Fst[5,]) #MAWR WIWR no Z 0.7544326
```

# Scans for introgressed loci
To scan for loci, I produced genotype-by-individual plots for each chromosome. These plots display the genotypes of each sample at high-Fst loci across a chromosome, displaying whether the sample carries the Pacific-associated allele, the Winter-associated allele, or is heterozygous. Hybrids are heterozygous at most of these loci, while un-introgressed samples are not heterozygous.

```R
# Genotype by individual plot ----
# Graph genotypes of high-Fst loci along a chromosome
# choose only loci that are variable in the dataset (SNPs), and above an Fst threshhold
groups.to.compare <- "PAWR_WIWR"   # can choose a pair or Fst_among for multi-group Fst_among

WC84_Fst.cutoff <- 0.9 #select a cutoff to display only SNPs above the cut-off
missing.fraction.max <- 0.4  # only show SNPs with less than this fraction of missing data among individuals
start.pos <- position.min
end.pos <- position.max
plot.groups <- groups
plot.group.colors <- group.colors
option <- 2    # choose whether to use color to show (1) ref/alt alleles or (2) group1/group2 alleles
# if option 2, then identify those with larger alternate allele freqs in group 1, and reverse the direction of genotypes (so ref and alt get switched--keep in mind for below):
group1 <- "PAWR"   # these groups will determine the color used in the graph
group2 <- "WIWR"
plotGenotypeByIndividual(groups.to.compare, WC84_Fst.cutoff, missing.fraction.max,
                         start.pos, end.pos, pos, WC84_Fst, combo.NApass, 
                         num_loc_cols, freqs, plot.groups, plot.group.colors,
                         option, group1, group2, chr)
```

Next, I was just interested to know how many SNPs total have Fst greater than 0.9 on a particular chromosome (whichever chromosome has had its data loaded into R)
```R
#How many SNPs are above Fst of 0.9 on chrZ?
#The values are kept in the variable WC84_Fst
groups.to.compare <- "PAWR_WIWR"  
WC84_Fst.cutoff <- 0.9 #select a cutoff to display only SNPs above the cut-off
row.choice <- which(rownames(WC84_Fst) == groups.to.compare)
selection <- (WC84_Fst[row.choice,] > WC84_Fst.cutoff) & (!is.na(WC84_Fst[row.choice,])) & (pos$position > start.pos) & (pos$position < end.pos) #can define start.pos and end.pos to select a smaller region of the chromosome
sum(selection) #this is just for one chromosome. 
#result: 299 SNPs on the Z.
```

# Fst per chromosome

Next, we can calculate the average Fst of each chromosome using the Weir + Cockerham 1984 method for multiple alleles.

```R

#get FST for entire chromosomes using Weir + Cockerham 1984 method for multiple alleles
chromosomes.to.analyze <- c("1", "1A", "2", "3","4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15","17","18","19","20","21","22","23","24","25","26","27","28","Z")

Analysis_set <- 1  # 1: all samples, only SNPs;    2: all samples, with invariant sites

if (Analysis_set == 1) {
  # choose path and filename for the 012NA files
  base.name <- "PAWR_WIWR_012NA_files/PAWR_WIWR.genotypes.SNPs_only"
  filename.text.middle <- ".max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 30 
  filename.text.end <-".MQ20.lowHet.tab"
  tag.name <- ".all_samples_abcfst."   # choose a tag name for this analysis
  # indicate name of metadata file, a text file with these column headings:
  # ID	location	group	Fst_group	plot_order
  metadata.file <- "PAWRWIWR_Metadata.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 77  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 1000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  #groups <- c("PAWR", "WIWR")
  #groups <- c("PAWR", "MAWR")
  groups <- c("WIWR", "MAWR")
  group.colors <- c("blue", "red")
  group_count <- length(groups)
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("PAWR_MAWR")
  group.colors.WC84_Fst <- c("purple")
} 

# Option to focus on a region of chromosome ----
focus.region <- F  # choose T for a subset of the chromosome, F for the whole thing)

calculate_or_load_stats <- 1  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) load per-site and windowed data from file

# Option to save the per-site stats
saveSiteInfo <- F # TRUE     # TRUE   # If TRUE, will save a file for per-site stats
saveWindowedStats <- F # TRUE     #TRUE   # If TRUE, will save a file for per-window stats
load.rolling.means <- F   #FALSE     #   FALSE    # If TRUE, will load rolling mean data (rather than calculate)
locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])
meanFst <- NA
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
  # combine metadata with genotype data
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  # If need to filter out individuals, based on low read number:
  filter <- F
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20,-163),]
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
    
    # calculate rownames for pairwise comparisons, for use in Dxy and Fst matrices:
    # rownames <- getPairwiseNames(groups)   # NOT NEEDED HERE since called from getDxy

    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    temp.list <- getWC84abc(freqs, sample_size, groups, among=FALSE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
    WC84_abc <- temp.list$WC84_abc
    aWC84 <- temp.list$aWC84
    bWC84 <- temp.list$bWC84
    cWC84 <- temp.list$cWC84
    rm(temp.list)
    print("Calculated WC84_abc values")
    
    if (saveSiteInfo == TRUE) {  # save the per-site stats, if chosen to in Intro section
      #save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R")) #uncorrected pi
      save(pos, freqs, sample_size, site_pi_nb, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
      print("Saved summary site stats")
      WC84_Fst_among <- WC84_Fst[rownames(WC84_Fst)=="Fst_among",]
      print(paste0(length(WC84_Fst_among)," markers in total (",sum(!is.na(WC84_Fst_among))," variable and ",sum(is.na(WC84_Fst_among)), " invariant)"))
    } else print("Site stats not saved")
    
  } else if (calculate_or_load_stats==2 | calculate_or_load_stats==3) {
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print("Loaded saved summary stats")
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
  
  #Here is the WC84 calculation of Fst as per pg 1364 of WC84 which works for multiple alleles
  aWC84[is.infinite(aWC84)] <- NA
  WC84_abc[is.infinite(WC84_abc)] <- NA
  chr_Fst <- (2*rowSums(aWC84, na.rm = T, dims = 1))/(2*rowSums(WC84_abc, na.rm = T, dims = 1))
  print(paste("Fst for ", chr, ":", chr_Fst))
  meanFst[i] <- chr_Fst
}   
# End main loop ----
#get Fst results
meanFst
```





















```R

# Option to focus on a region of chromosome ----
focus.region <- T  # choose T for a subset of the chromosome, F for the whole thing)
if (focus.region==T){
  position.min <- 0000000
  position.max <- 3000000
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
  filter <- F
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20,-163),]
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
    #site_pi <- getSitePi(freqs) #for uncorrected pi, biased low at low sample size
    site_pi_nb <- getSitePi_nb(freqs, sample_size) #corrected pi
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
      #save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R")) #uncorrected pi
      save(pos, freqs, sample_size, site_pi_nb, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
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
    #temp.list <- getWindowedPi(site_pi, pos, window_size, step_size) #uncorrected pi
    temp.list <- getWindowedPi(site_pi_nb, pos, window_size, step_size)
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


###This section performs another PCA plot for pacificus vs WIWR PCA#####
## Make PCA plot ----
# choose only loci that are variable in the dataset (SNPs), and (optionally) above an Fst threshhold
# groups and colors determined in Intro section, under groups.to.plot.PCA and group.colors.PCA
Fst.filter <- F
Fst.cutoff <- 0.01
# choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Fst_among"
#groups.to.compare <- "nigrifrons_auduboni"
axes <- 3

groups.to.plot.PCA <- c("PAWR", "WIWR", "Hybrid")
group.colors.PCA <- c("blue", "red", "purple")

PCA_results <- plotPCA(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                       groups.to.plot.PCA, group.colors.PCA, axes, flip1=F, flip2=F)

```

It was suggested to check whether fst is correlated with window size in this dataset. 
```R
getWindowSize <- function(pos, window_size, step_size){
  rolling.min.pos <- rollapply(pos$position, width=window_size, FUN=min, by=step_size, align="left")
  rolling.max.pos <- rollapply(pos$position, width=window_size, FUN=max, by=step_size, align="left")
  rolling.size <- rolling.max.pos - rolling.min.pos
  return(rolling.size)
}
rolling.size <- getWindowSize(pos, window_size, step_size)

mean_WC84_Fst=rolling.mean.WC84_Fst[1,]
mean_Dxy=rolling.mean.Dxy[1,]
mean_pi=rolling.mean.pi[1,]
ggplot(as.data.frame(mean_WC84_Fst), aes(x=rolling.size, y=mean_WC84_Fst))+
  geom_point()+
  theme_classic()+
  xlab("Window Size")+
  ylab("Fst")+
  stat_smooth(method=lm, aes(x=rolling.size, y=mean_WC84_Fst))
ggplot(as.data.frame(mean_pi), aes(x=rolling.size, y=mean_pi))+
  geom_point()+
  theme_classic()+
  xlab("Window Size")+
  ylab("pi_w")+
  stat_smooth(method=lm, aes(x=rolling.size, y=mean_pi))
ggplot(as.data.frame(mean_Dxy), aes(x=rolling.size, y=mean_Dxy))+
  geom_point()+
  theme_classic()+
  xlab("Window Size")+
  ylab("pi_b")+
  stat_smooth(method=lm, aes(x=rolling.size, y=mean_Dxy))



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
    combo.NApass.whole.chr <- combo[c(-20,-163),]
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

  load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
  print("Loaded saved summary stats")
  
  load(paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
   
  # Make plot for chr ----
  # make ggplots for quick inspection of rolling mean results
  rolling.size <- getWindowSize(pos, window_size, step_size)

  mean_WC84_Fst=rolling.mean.WC84_Fst[1,]
  mean_Dxy=rolling.mean.Dxy[1,]
  mean_pi=rolling.mean.pi[1,]
  ggplot(as.data.frame(mean_WC84_Fst), aes(x=rolling.size, y=mean_WC84_Fst))+
    geom_point()+
    theme_classic()+
    xlab("Window Size")+
    ylab("Fst")+
    stat_smooth(method=lm, aes(x=rolling.size, y=mean_WC84_Fst))
  ggplot(as.data.frame(mean_pi), aes(x=rolling.size, y=mean_pi))+
    geom_point()+
    theme_classic()+
    xlab("Window Size")+
    ylab("pi_w")+
    stat_smooth(method=lm, aes(x=rolling.size, y=mean_pi))
  ggplot(as.data.frame(mean_Dxy), aes(x=rolling.size, y=mean_Dxy))+
    geom_point()+
    theme_classic()+
    xlab("Window Size")+
    ylab("pi_b")+
    stat_smooth(method=lm, aes(x=rolling.size, y=mean_Dxy))

} 
```
