# Introduction ----
# genomics_R_functions.R
# Written by Darren Irwin in 2015-2016, initially for the purpose of analyzing GBS data
# from Greenish Warblers.
# This has been modified a bit for use in North American warbler analyses.
# These functions are used in the analysis presented in the following paper:
#
# Irwin DE, Alcaide M, Delmore KE, Irwin JH, Owens GL. Recurrent selection explains parallel evolution of genomic regions of high relative but low absolute differentiation in a ring species. In revision, Molecular Ecology.
#
# also posted on bioRxiv:
#
# Irwin DE, Alcaide M, Delmore KE, Irwin JH, Owens GL. 2016. Recurrent selection explains parallel evolution of genomic regions of high relative but low absolute differentiation in greenish warblers. bioRxiv, doi: http://dx.doi.org/10.1101/041467
#
# To run the actual analyses in the above paper, use the accompanying script file:
# # GW_GBS_R_analysis_script_for_Dryad.R
#
# If you use these scripts / functions, please cite the above paper and/or the Dryad package.

# Libraries ----
# Need to first install these by typing install.packages("[[insert library name]]")
library(zoo)
library(grid)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(fields)
# For pcaMethods, need to do these commands:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("pcaMethods")
library(pcaMethods)  # see note 4 lines up to install this

# Functions ----

# Calculate allele freqs and sample sizes (use column Fst_group)
getFreqsAndSampleSizes <- function(dataset, num_intro_cols, groups){
  group_count <- length(groups)
  freqs <- matrix(nrow = group_count, ncol = (length(dataset[1,]) - num_intro_cols), dimnames = list(groups, colnames(dataset)[(num_intro_cols+1):length(dataset[1,])]))
  sample_size <- matrix(nrow = group_count, ncol = (length(dataset[1,]) - num_intro_cols), dimnames = list(groups, colnames(dataset)[(num_intro_cols+1):length(dataset[1,])]))
  for (i in 1:group_count) {
    selection <- dataset[(dataset$Fst_group == groups[i]),]
    sums <- colSums(selection[,(num_intro_cols+1):length(selection[1,])], na.rm=TRUE)
    counts <- colSums(!is.na(selection[,(num_intro_cols+1):length(selection[1,])]))  #counts rows per column that are not NA
    freqs[i,] <- 0.5 * sums / counts
    sample_size[i,] <- counts
  }
  return(list(freqs=freqs, sample_size=sample_size))
}


# calculate nucleotide diversity (pi) at each site for each population,
# using allele freqs as input.
# (note this is biased low at limited sample size)
getSitePi <- function(freqs){
  site_pi <- 2*freqs*(1-freqs)
  return(site_pi)
}

# calculate non-biased nucleotide diversity (pi) at each site for each population,
# with correction for sample size,
# using allele freqs and sample size as input:
getSitePi_nb <- function(freqs, sample_size){
  site_pi_nb <- (2*sample_size/(2*sample_size-1)) * 2*freqs*(1-freqs)
  # change infinites (when sample_size is 1) to NaN:
  site_pi_nb[is.infinite(site_pi_nb)] <- NA
  return(site_pi_nb) 
}

# get row names for use in Dxy and Fst matrices, 
# or any others that compare populations pairwise:
getPairwiseNames <- function(groups){
  group_count <- length(groups)
  pairwise_names <- vector(mode="character", length = choose(group_count, 2))
  rowcount <- 1
  for (i in 1:(group_count-1)) {
    for (j in (i+1):group_count) {
      pairwise_names[rowcount] <- paste(groups[i], groups[j], sep="_")
      rowcount <- rowcount + 1
    }	
  }
  return(pairwise_names)
}


# calculate pairwise Dxy per site, using data in "freqs" and groups in "groups"
getDxy <- function(freqs, groups){
  rownames <- getPairwiseNames(groups)
  group_count <- length(groups)
  Dxy <- matrix(nrow = length(rownames), ncol = length(freqs[1,]), dimnames = list(rownames, colnames(freqs)))
  rowcount <- 1
  for (i in 1:(group_count-1)) {
    for (j in (i+1):group_count) {
      Dxy[rowcount,] <- freqs[i,]*(1-freqs[j,]) + freqs[j,]*(1-freqs[i,])
      rowcount <- rowcount + 1
    }	# the above "j" loop calculates Dxy per site
  }
  return(Dxy)
}


# calculate Fst (and numerator and denominator) for each site (bp), 
# between pairs of groups (so pops (r) is 2), 
# using the Weir&Cockerham 1984 approach to correct for sample size and number of pops.
# Set "among=TRUE" if wanting it to include a line for Fst among all groups.
getWC84Fst <- function(freqs, sample_size, groups, among=FALSE){
  rownames <- getPairwiseNames(groups)
  group_count <- length(groups)
  WC84_Fst <- matrix(nrow = length(rownames), ncol = length(freqs[1,]), dimnames = list(rownames, colnames(freqs)))
  WC84_Fst_numerator <- matrix(nrow = length(rownames), ncol = length(freqs[1,]), dimnames = list(rownames, colnames(freqs)))
  WC84_Fst_denominator <- matrix(nrow = length(rownames), ncol = length(freqs[1,]), dimnames = list(rownames, colnames(freqs)))
  rowcount <- 1
  for (i in 1:(group_count-1)) {
    for (j in (i+1):group_count) {
      mean_sample_size <- (sample_size[i,] + sample_size[j,]) / 2
      mean_freqs <- (sample_size[i,]*freqs[i,] + sample_size[j,]*freqs[j,]) / (2*mean_sample_size)
      CV_sample_size <- sqrt(((sample_size[i,]-mean_sample_size)^2 + (sample_size[j,]-mean_sample_size)^2)/1) / mean_sample_size     # Note that used the equation for sample SD here--Mike & I solved the n_c equation in Weir&Cockerham1984 and confirmed this. (see the division by 1 in the numerator, prior to square root)
      allele_freq_variance <- (sample_size[i,]*((freqs[i,]-mean_freqs)^2) + sample_size[j,]*((freqs[j,]-mean_freqs)^2)) / ((2-1)*mean_sample_size)   # using sample variance, as per W&C84
      WC84_Fst_numerator[rowcount,] <- allele_freq_variance - ((mean_freqs*(1-mean_freqs) - (1/2)*allele_freq_variance)/(2*mean_sample_size - 1))
      WC84_Fst_denominator.term1 <- (1 - (2*mean_sample_size*((CV_sample_size)^2)/((2*mean_sample_size - 1)*2))) * mean_freqs*(1-mean_freqs) 
      WC84_Fst_denominator.term2 <- (1 + ((2*mean_sample_size*(CV_sample_size)^2)/((2*mean_sample_size-1)*2))) * allele_freq_variance / 2
      WC84_Fst_denominator[rowcount,] <- WC84_Fst_denominator.term1 + WC84_Fst_denominator.term2
      WC84_Fst[rowcount,] <- WC84_Fst_numerator[rowcount,] / WC84_Fst_denominator[rowcount,]		
      rowcount <- rowcount + 1
    }
  }
  if (among==TRUE){
    # Add a line for the "Fst_among" for each site, 
    # among all the groups listed in "groups" variable above (so r can be greater than 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size.
    WC84_Fst_among <- matrix(nrow = 1, ncol = length(freqs[1,]), dimnames = list("Fst_among", colnames(freqs)))
    WC84_Fst_among_numerator <- matrix(nrow = 1, ncol = length(freqs[1,]), dimnames = list("Fst_among", colnames(freqs)))
    WC84_Fst_among_denominator <- matrix(nrow = 1, ncol = length(freqs[1,]), dimnames = list("Fst_among", colnames(freqs)))
    for (i in 1:length(freqs[1,])) {
      pop.freqs <- freqs[,i]
      pop.freqs.nan_removed <- pop.freqs[is.finite(pop.freqs)]
      num.pops <- length(pop.freqs.nan_removed)    # "r" in Weir&Cockerham
      pop.sample_sizes <- sample_size[,i]
      pop.sample_sizes.positive <- pop.sample_sizes[is.finite(pop.freqs)]
      mean_sample_size <- mean(pop.sample_sizes.positive)
      mean_freq <- mean(pop.freqs.nan_removed)
      CV_sample_size <- sqrt(sum((pop.sample_sizes.positive - mean_sample_size)^2) / (num.pops - 1)) / mean_sample_size     # Note that used the equation for sample SD here--Mike & I solved the n_c equation in Weir&Cockerham1984 and confirmed this.
      allele_freq_variance <- sum(pop.sample_sizes.positive*((pop.freqs.nan_removed - mean_freq)^2)) / ((num.pops - 1)*mean_sample_size)    # using sample variance, as per W&C84
      WC84_Fst_among_numerator[1,i] <- allele_freq_variance - ((mean_freq*(1-mean_freq) - ((num.pops-1)/num.pops)*allele_freq_variance)/(2*mean_sample_size - 1))
      WC84_Fst_among_denominator.term1 <- (1 - (2*mean_sample_size*((CV_sample_size)^2)/((2*mean_sample_size - 1)*num.pops))) * mean_freq*(1-mean_freq) 
      WC84_Fst_among_denominator.term2 <- (1 + ((2*mean_sample_size*(num.pops-1)*(CV_sample_size^2))/((2*mean_sample_size-1)*num.pops))) * allele_freq_variance / num.pops
      WC84_Fst_among_denominator[1,i] <- WC84_Fst_among_denominator.term1 + WC84_Fst_among_denominator.term2
      WC84_Fst_among[1,i] <- WC84_Fst_among_numerator[1,i] / WC84_Fst_among_denominator[1,i]		
    }
    WC84_Fst <- rbind(WC84_Fst, WC84_Fst_among)
    WC84_Fst_numerator <- rbind(WC84_Fst_numerator, WC84_Fst_among_numerator)
    WC84_Fst_denominator <- rbind(WC84_Fst_denominator, WC84_Fst_among_denominator)
  }
  return(list(WC84_Fst=WC84_Fst, WC84_Fst_numerator=WC84_Fst_numerator, WC84_Fst_denominator=WC84_Fst_denominator))
}


# calculate windowed pi, in whole windows starting on left side of chromosome
getWindowedPi <- function(site_pi, pos, window_size, step_size){
  rolling.mean.pos.pi <- rollapply(pos$position, width=window_size, FUN=mean, by=step_size, align="left")
  rolling.mean.pi <- matrix(nrow = length(groups), ncol = length(rolling.mean.pos.pi), dimnames = list(rownames(site_pi), round(rolling.mean.pos.pi)))
  for (i in 1:length(site_pi[,1])) {
    rolling.mean.pi[i,] <- rollapply(site_pi[i,], width=window_size, FUN=mean, na.rm = TRUE, by=step_size, align="left")
  } 
  return(list(rolling.mean.pos.pi=rolling.mean.pos.pi, rolling.mean.pi=rolling.mean.pi))
}


# calculate windowed Dxy, in whole windows starting on left side of chromosome
getWindowedDxy <- function(Dxy, pos, window_size, step_size){
  rolling.mean.pos.Dxy <- rollapply(pos$position, width=window_size, FUN=mean, by=step_size, align="left")
  rolling.mean.Dxy <- matrix(nrow = length(Dxy[,1]), ncol = length(rolling.mean.pos.Dxy), dimnames = list(rownames(Dxy), round(rolling.mean.pos.Dxy)))
  for (i in 1:length(Dxy[,1])) {
    rolling.mean.Dxy[i,] <- rollapply(Dxy[i,], width=window_size, FUN=mean, na.rm = TRUE, by=step_size, align="left")
  }
  return(list(rolling.mean.pos.Dxy=rolling.mean.pos.Dxy, rolling.mean.Dxy=rolling.mean.Dxy))
}


# calculate windowed Fst according to according to Weir&Cockerham1984 (with sample size and pop number correction),
# calculated as windowed numerator over windowed denominator,
# in whole windows starting on left side of chromosome
getWindowedWC84_Fst <- function(WC84_Fst_numerator, WC84_Fst_denominator, pos, window_size, step_size){
  rolling.mean.pos.WC84_Fst <- rollapply(pos$position, width=window_size, FUN=mean, by=step_size, align="left")   # if want partial window on right, add: partial=window_size/2
  rolling.mean.WC84_Fst <- matrix(nrow = length(WC84_Fst_numerator[,1]), ncol = length(rolling.mean.pos.WC84_Fst), dimnames = list(rownames(WC84_Fst_numerator), round(rolling.mean.pos.WC84_Fst)))
  for (i in 1:length(WC84_Fst_numerator[,1])) {
    rolling.mean.WC84_Fst_numerator <- rollapply(WC84_Fst_numerator[i,], width=window_size, FUN=mean, na.rm = TRUE, by=step_size, align="left")
    rolling.mean.WC84_Fst_denominator <- rollapply(WC84_Fst_denominator[i,], width=window_size, FUN=mean, na.rm = TRUE, by=step_size, align="left")
    rolling.mean.WC84_Fst[i,] <- rolling.mean.WC84_Fst_numerator / rolling.mean.WC84_Fst_denominator
  }
  return(list(rolling.mean.pos.WC84_Fst=rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst=rolling.mean.WC84_Fst))
}


# make plots for quick inspection of rolling mean results
makeRollingMeanPlots <- function(rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst, 
                                 rolling.mean.pos.Dxy, rolling.mean.Dxy, 
                                 rolling.mean.pos.pi, rolling.mean.pi,
                                 group.colors.pi, groups.to.plot.pi,
                                 group.colors.Dxy, groups.to.plot.Dxy,
                                 group.colors.WC84_Fst, groups.to.plot.WC84_Fst,
                                 region.text){
  # arrange data for ggplot:
  mean.WC84_Fst.data <- data.frame(mean_position=rolling.mean.pos.WC84_Fst, mean_WC84_Fst=rolling.mean.WC84_Fst[1,], populations=rep(rownames(rolling.mean.WC84_Fst)[1], length(rolling.mean.pos.WC84_Fst)))
  if (length(rolling.mean.Dxy[,1]) > 1) {
    for (i in 2:length(rolling.mean.WC84_Fst[,1])) {
    mean.WC84_Fst.data <- rbind(mean.WC84_Fst.data, data.frame(mean_position=rolling.mean.pos.WC84_Fst, mean_WC84_Fst=rolling.mean.WC84_Fst[i,], populations=rep(rownames(rolling.mean.WC84_Fst)[i], length(rolling.mean.pos.WC84_Fst))))
    }
  }
  mean.Dxy.data <- data.frame(mean_position=rolling.mean.pos.Dxy, mean_Dxy=rolling.mean.Dxy[1,], populations=rep(rownames(rolling.mean.Dxy)[1], length(rolling.mean.pos.Dxy)))
  if (length(rolling.mean.Dxy[,1]) > 1) {
    for (i in 2:length(rolling.mean.Dxy[,1])) {
      mean.Dxy.data <- rbind(mean.Dxy.data, data.frame(mean_position=rolling.mean.pos.Dxy, mean_Dxy=rolling.mean.Dxy[i,], populations=rep(rownames(rolling.mean.Dxy)[i], length(rolling.mean.pos.Dxy))))
    }
  }
  mean.pi.data <- data.frame(mean_position=rolling.mean.pos.pi, mean_pi=rolling.mean.pi[1,], taxon=rep(rownames(rolling.mean.pi)[1], length(rolling.mean.pos.pi)))
  for (i in 2:length(rolling.mean.pi[,1])) {
    mean.pi.data <- rbind(mean.pi.data, data.frame(mean_position=rolling.mean.pos.pi, mean_pi=rolling.mean.pi[i,], taxon=rep(rownames(rolling.mean.pi)[i], length(rolling.mean.pos.pi))))
  }
  # define ggplots:
  plot_data <- mean.pi.data[mean.pi.data$taxon %in% groups.to.plot.pi,]
  plot_data$taxon <- factor(plot_data$taxon, levels=groups.to.plot.pi)
  pi_plot <- ggplot(data=plot_data, aes(x=mean_position, y=mean_pi, colour=taxon)) + 
    geom_line() + 
    scale_colour_manual(values = group.colors.pi) +
    theme(legend.key.size = unit(0.4, "cm"), legend.text=element_text(size=8))	
  
  plot_data <- mean.Dxy.data[mean.Dxy.data$populations %in% groups.to.plot.Dxy,]
  plot_data$populations <- factor(plot_data$populations, levels=groups.to.plot.Dxy)
  Dxy_plot <- ggplot(data=plot_data, aes(x=mean_position, y=mean_Dxy, colour=populations)) +
    geom_line() +
    scale_colour_manual(values = group.colors.Dxy) +
    theme(legend.key.size = unit(0.4, "cm"), legend.text=element_text(size=8))
  
  plot_data <- mean.WC84_Fst.data[mean.WC84_Fst.data$populations %in% groups.to.plot.WC84_Fst,]
  plot_data$populations <- factor(plot_data$populations, levels=groups.to.plot.WC84_Fst)
  WC84_Fst_plot <- ggplot(data=plot_data, aes(x=mean_position, y=mean_WC84_Fst, colour=populations)) +
    geom_line() +
    scale_colour_manual(values = group.colors.WC84_Fst) +
    theme(legend.key.size = unit(0.4, "cm"), legend.text=element_text(size=8))
  
  # Generate the plot figure:
  quartz(title=paste0(region.text,": pi, Fst, and Dxy"), width=10, height=6)
  gA <- ggplot_gtable(ggplot_build(WC84_Fst_plot))
  gB <- ggplot_gtable(ggplot_build(Dxy_plot))    # Dxy_plot
  gC <- ggplot_gtable(ggplot_build(pi_plot))		# pi_plot
  gA$widths <- gB$widths
  gC$widths <- gB$widths
  grid.newpage()
  grid.arrange(gA, gB, gC, nrow = 3)
}


# Pi graph: show both site data and rolling means
# This works best with three groups
makePiPlot <- function(pos, site_pi, region.text,
                       rolling.mean.pos.pi, rolling.mean.pi,
                       groups.to.plot.pi, group.colors.pi){
  # For group 1
  graph.data <- data.frame(position=pos$position, pi=site_pi[rownames(site_pi)==groups.to.plot.pi[1],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.pi, pi=amp*rolling.mean.pi[rownames(rolling.mean.pi)==groups.to.plot.pi[1],], point_type=rep("windows", length(rolling.mean.pos.pi)))) 
  plot.pi.group1 <- ggplot(graph.data, aes(x=position, y=pi, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype=point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.pi[1])) + scale_linetype_manual(values=c("blank", "solid"))
  # For group 2
  graph.data <- data.frame(position=pos$position, pi=site_pi[rownames(site_pi)==groups.to.plot.pi[2],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.pi, pi=amp*rolling.mean.pi[rownames(rolling.mean.pi)==groups.to.plot.pi[2],], point_type=rep("windows", length(rolling.mean.pos.pi)))) 
  plot.pi.group2 <- ggplot(graph.data, aes(x=position, y=pi, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype=point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.pi[2])) + scale_linetype_manual(values=c("blank", "solid"))
  # For group 3
  graph.data <- data.frame(position=pos$position, pi=site_pi[rownames(site_pi)==groups.to.plot.pi[3],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.pi, pi=amp*rolling.mean.pi[rownames(rolling.mean.pi)==groups.to.plot.pi[3],], point_type=rep("windows", length(rolling.mean.pos.pi)))) 
  plot.pi.group3 <- ggplot(graph.data, aes(x=position, y=pi, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype=point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.pi[3])) + scale_linetype_manual(values=c("blank", "solid"))
  # Generate the plot figure for pi:
  string <- paste(groups.to.plot.pi, collapse=", ")
  quartz(title=paste(region.text,": pi in ", string, sep=""), width=10, height=5)
  gA <- ggplot_gtable(ggplot_build(plot.pi.group1))
  gB <- ggplot_gtable(ggplot_build(plot.pi.group2))
  gC <- ggplot_gtable(ggplot_build(plot.pi.group3))
  gA$widths <- gB$widths
  gC$widths <- gB$widths
  grid.newpage()
  grid.arrange(gA, gB, gC, nrow = 3) 
}


# Dxy graph: show both site data and rolling means
# This works best with three groups
makeDxyPlot <- function(pos, Dxy, region.text, 
                        rolling.mean.pos.Dxy, rolling.mean.Dxy,
                        groups.to.plot.Dxy, group.colors.Dxy){
  # For group 1
  graph.data <- data.frame(position=pos$position, Dxy=Dxy[rownames(Dxy)==groups.to.plot.Dxy[1],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.Dxy, Dxy= amp*rolling.mean.Dxy[rownames(rolling.mean.Dxy)==groups.to.plot.Dxy[1],], point_type=rep("windows", length(rolling.mean.pos.Dxy))))
  plot.Dxy.group1 <- ggplot(graph.data, aes(x=position, y=Dxy, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype= point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.Dxy[1])) + scale_linetype_manual(values=c("blank", "solid"))
  # For group 2
  graph.data <- data.frame(position=pos$position, Dxy=Dxy[rownames(Dxy)==groups.to.plot.Dxy[2],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.Dxy, Dxy= amp*rolling.mean.Dxy[rownames(rolling.mean.Dxy)==groups.to.plot.Dxy[2],], point_type=rep("windows", length(rolling.mean.pos.Dxy))))
  plot.Dxy.group2 <- ggplot(graph.data, aes(x=position, y=Dxy, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype= point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.Dxy[2])) + scale_linetype_manual(values=c("blank", "solid"))
  # For group 3
  graph.data <- data.frame(position=pos$position, Dxy=Dxy[rownames(Dxy)==groups.to.plot.Dxy[3],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.Dxy, Dxy= amp*rolling.mean.Dxy[rownames(rolling.mean.Dxy)==groups.to.plot.Dxy[3],], point_type=rep("windows", length(rolling.mean.pos.Dxy))))
  plot.Dxy.group3 <- ggplot(graph.data, aes(x=position, y=Dxy, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype= point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.Dxy[3])) + scale_linetype_manual(values=c("blank", "solid"))
  # Generate the plot figure for Dxy:
  string <- paste(groups.to.plot.Dxy, collapse=", ")
  quartz(title=paste0(region.text,": Dxy between ", string, sep=""), width=10, height=5)
  gA <- ggplot_gtable(ggplot_build(plot.Dxy.group1))
  gB <- ggplot_gtable(ggplot_build(plot.Dxy.group2))
  gC <- ggplot_gtable(ggplot_build(plot.Dxy.group3))
  gA$widths <- gB$widths
  gC$widths <- gB$widths
  grid.newpage()
  grid.arrange(gA, gB, gC, nrow = 3)
}


# Fst graph: show both site data and rolling means
# Currently uses WC84_Fst but would work with others
# This works best with three groups
makeFstPlot <- function(pos, WC84_Fst, region.text, 
                        rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst,
                        groups.to.plot.WC84_Fst, group.colors.WC84_Fst){
  # For group 1
  graph.data <- data.frame(position=pos$position, WC84_Fst=WC84_Fst[rownames(WC84_Fst)==groups.to.plot.WC84_Fst[1],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.WC84_Fst, WC84_Fst= amp*rolling.mean.WC84_Fst[rownames(rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[1],], point_type=rep("windows", length(rolling.mean.pos.WC84_Fst))))
  plot.WC84_Fst.group1 <- ggplot(graph.data, aes(x=position, y=WC84_Fst, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype= point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.WC84_Fst[1])) + scale_linetype_manual(values=c("blank", "solid"))
  # For group 2
  graph.data <- data.frame(position=pos$position, WC84_Fst=WC84_Fst[rownames(WC84_Fst)==groups.to.plot.WC84_Fst[2],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.WC84_Fst, WC84_Fst= amp*rolling.mean.WC84_Fst[rownames(rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[2],], point_type=rep("windows", length(rolling.mean.pos.WC84_Fst))))
  plot.WC84_Fst.group2 <- ggplot(graph.data, aes(x=position, y=WC84_Fst, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype= point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.WC84_Fst[2])) + scale_linetype_manual(values=c("blank", "solid"))
  # For group 3
  graph.data <- data.frame(position=pos$position, WC84_Fst=WC84_Fst[rownames(WC84_Fst)==groups.to.plot.WC84_Fst[3],], point_type=rep("sites", length(pos$position)))
  graph.data <- rbind(graph.data, data.frame(position=rolling.mean.pos.WC84_Fst, WC84_Fst= amp*rolling.mean.WC84_Fst[rownames(rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[3],], point_type=rep("windows", length(rolling.mean.pos.WC84_Fst))))
  plot.WC84_Fst.group3 <- ggplot(graph.data, aes(x=position, y=WC84_Fst, group=point_type)) + geom_point(aes(color=point_type)) + geom_line(aes(linetype= point_type, color=point_type)) + scale_color_manual(values=c("black",group.colors.WC84_Fst[3])) + scale_linetype_manual(values=c("blank", "solid"))
  
  # Generate the plot figure for WC84_Fst:
  string <- paste(groups.to.plot.WC84_Fst, collapse=", ")
  quartz(title=paste0(region.text,": WC84_Fst between ", string, sep=""), width=10, height=5)
  gA <- ggplot_gtable(ggplot_build(plot.WC84_Fst.group1))
  gB <- ggplot_gtable(ggplot_build(plot.WC84_Fst.group2))
  gC <- ggplot_gtable(ggplot_build(plot.WC84_Fst.group3))
  gA$widths <- gB$widths
  gC$widths <- gB$widths
  grid.newpage()
  grid.arrange(gA, gB, gC, nrow = 3)
}


# Detailed chr plot
# Plots of Fst, Dxy, and pi for three groups across chromosome
# Note there is a small amount of vertical jitter in point locations
makeDetailedChrPlots <- function(pos, 
                                 WC84_Fst, rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst,
                                 Dxy, rolling.mean.pos.Dxy, rolling.mean.Dxy,
                                 site_pi, rolling.mean.pos.pi, rolling.mean.pi, 
                                 group.colors.pi, groups.to.plot.pi,
                                 group.colors.Dxy, groups.to.plot.Dxy,
                                 group.colors.WC84_Fst, groups.to.plot.WC84_Fst,
                                 region.text) {
  window.width <- 8
  window.height <- 10
  # numbers below are in inches:
  left.margin <- 1 / window.width  # the left outer margin of the figure (number is inches)
  right.margin <- 1.5 / window.width
  top.margin <- 0.5 / window.height
  bottom.margin <- 0.5 / window.height
  row.height <- 2.5 / window.height  # height of each row of the figure (Fst row, Dxy row, and pi row)
  plot.height.Fst <- row.height / length(groups.to.plot.WC84_Fst)
  plot.height.Dxy <- row.height / length(groups.to.plot.Dxy)
  plot.height.pi <- row.height / length(groups.to.plot.pi)
  row.space <- 0.25 / window.height  # gap between rows 
  #bp.per.inch <- 20000000  # number of bp of sequence per inch
  data.bp.length <- pos$position[length(pos$position)]
  # start.bp <- 0
  # end.bp <- data.bp.length
  start.bp <- pos$position[1]      #53500000
  end.bp <- pos$position[length(pos$position)]      #55500000
  fig.left.loc <- left.margin  
  fig.right.loc <- 1 - right.margin
  
  amp.Fst <- 1
  amp.Dxy <- 10
  amp.pi <- 10
  linewidth <- 4
  point.char <- 16
  pointsize <- 0.25
  line.transparency <- 0.8
  jitter.Fst <- 0.025
  jitter.Dxy <- 0.025
  jitter.pi <- 0.0125
  
  quartz(title=paste0("WC84_Fst, Dxy, and pi across ", region.text), width=window.width, height=window.height)  # width and height in inches
  label.size <- 1
  new.fig.parameter <- FALSE
  row.num <- 1
  for (j in 1:length(groups.to.plot.WC84_Fst)) {
    par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+j*plot.height.Fst+(row.num-1)*(row.height+row.space)), 1-(top.margin+(j-1)*plot.height.Fst+(row.num-1)*(row.height+row.space))), new = new.fig.parameter, mai=c(0,0,0,0), cex=0.5)  # mai sets the margins of the plot to zero, matching
    new.fig.parameter <- TRUE
    plot(x=NULL, y=NULL, xlim=c(start.bp, end.bp), ylim=c(0, 1.2), yaxp=c(0,1, n=1), mgp=c(2,0.5,0), ylab=NA, xlab=NA, xaxt='n', yaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex=1, cex.lab=1.25, cex.axis=0.75)
    mtext(side=4, text=expression("Average " * italic('F')[ST]), col=group.colors.WC84_Fst[j], cex=1, las=1, line=1.5)
    mtext(side=2, text=expression(italic('F')[ST]), col="black", cex=1.2, las=1, line=2)
    axis(side=2, at=c(0,1), tcl=-0.25, las=1, col="black", col.ticks="black", cex=0.75, mgp=c(2,0.5,0))
    axis(side=4, at=c(0,1), tcl=-0.25, las=1, col="black", col.ticks=group.colors.WC84_Fst[j], cex=0.75, mgp=c(2,0.5,0))
    # for graph only, convert all negative Fst's to zero:
    Fst_for_graph <- WC84_Fst[rownames(WC84_Fst)==groups.to.plot.WC84_Fst[j],]
    Fst_for_graph[Fst_for_graph < 0] <- 0
    points(pos$position, jitter(Fst_for_graph, amount=jitter.Fst), pch=point.char, cex=pointsize)
    color.rgb <- col2rgb(group.colors.WC84_Fst[j])/255
    lines(rolling.mean.pos.WC84_Fst, amp.Fst*rolling.mean.WC84_Fst[rownames(rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[j],], col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line.transparency), lwd=linewidth)
  }
  # plot Dxy:
  row.num <- 2
  for (j in 1:length(groups.to.plot.Dxy)) {
    par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+j*plot.height.Dxy+(row.num-1)*(row.height+row.space)), 1-(top.margin+(j-1)*plot.height.Dxy+(row.num-1)*(row.height+row.space))), new = new.fig.parameter, mai=c(0,0,0,0), cex=0.5)  # mai sets the margins of the plot to zero, matching
    new.fig.parameter <- TRUE
    plot(x=NULL, y=NULL, xlim=c(start.bp, end.bp), ylim=c(0, 1.2), yaxp=c(0, 1, n=1), mgp=c(2,0.5,0), ylab=NA, xlab=NA, xaxt='n', yaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex=1, cex.lab=1.25, cex.axis=0.75)
    mtext(side=4, text=expression("10 * mean " * italic('D')[xy]), col=group.colors.Dxy[j], cex=1, las=1, line=1.5)
    mtext(side=2, text=expression(italic('D')[xy]), col="black", cex=1.2, las=1, line=2)
    axis(side=2, at=c(0,1), tcl=-0.25, las=1, col="black", col.ticks="black", cex=0.75, mgp=c(2,0.5,0))
    axis(side=4, at=c(0,1), tcl=-0.25, las=1, col="black", col.ticks=group.colors.Dxy[j], cex=0.75, mgp=c(2,0.5,0))
    Dxy_for_graph <- Dxy[rownames(Dxy)==groups.to.plot.Dxy[j],]
    points(pos$position, jitter(Dxy_for_graph, amount=jitter.Dxy), pch=point.char, cex=pointsize)
    color.rgb <- col2rgb(group.colors.Dxy[j])/255
    lines(rolling.mean.pos.Dxy, amp.Dxy*rolling.mean.Dxy[rownames(rolling.mean.Dxy)==groups.to.plot.Dxy[j],], col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line.transparency), lwd=linewidth)
  }
  # plot pi: 
  row.num <- 3
  for (j in 1:length(groups.to.plot.pi)) {
    par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+j*plot.height.pi+(row.num-1)*(row.height+row.space)), 1-(top.margin+(j-1)*plot.height.pi+(row.num-1)*(row.height+row.space))), new = new.fig.parameter, mai=c(0,0,0,0), cex=0.5)  # mai sets the margins of the plot to zero, matching
    new.fig.parameter <- TRUE
    plot(x=NULL, y=NULL, xlim=c(start.bp, end.bp), ylim=c(0, 0.7), yaxp=c(0, 0.5, n=1), mgp=c(2,0.5,0), ylab=NA, xlab=NA, xaxt='n', yaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex=1, cex.lab=1.25, cex.axis=0.75)
    mtext(side=4, text=expression("10 * mean " * italic(pi)), col=group.colors.pi[j], cex=1, las=1, line=1.5)
    mtext(side=2, text=expression(italic(pi) * " "), col="black", cex=1.2, las=1, line=2)
    axis(side=2, at=c(0,0.5), tcl=-0.25, las=1, col="black", col.ticks="black", cex=0.75, mgp=c(2,0.5,0))
    axis(side=4, at=c(0,0.5), tcl=-0.25, las=1, col="black", col.ticks=group.colors.pi[j], cex=0.75, mgp=c(2,0.5,0))
    site_pi_for_graph <- site_pi[rownames(site_pi)==groups.to.plot.pi[j],]
    points(pos$position, jitter(site_pi_for_graph, amount=jitter.pi), pch=point.char, cex=pointsize)
    color.rgb <- col2rgb(group.colors.pi[j])/255
    lines(rolling.mean.pos.pi, amp.pi*rolling.mean.pi[rownames(rolling.mean.pi)==groups.to.plot.pi[j],], col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line.transparency), lwd=linewidth)
  }
  #mtext(paste0("Location along chromosome ", region.text, sep=""), side=1, line=3)  #, adj=0, at=c(0.5,0.5))
  mtext(paste0("Location along ", region.text), side=1, line=3)  #, adj=0, at=c(0.5,0.5))
}


# Make PCA plot
# choose only loci that are variable in the dataset (SNPs), and (optionally) above an Fst threshhold
# groups and colors determined in Intro section, under groups.to.plot.PCA and group.colors.PCA
plotPCA <- function(Fst.filter, Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                    groups.to.plot.PCA, group.colors.PCA, axes, flip1, flip2) {
  #choose the comparison to apply the Fst filter to in the command below:
  if (Fst.filter == TRUE) {
    row.choice <- which(rownames(WC84_Fst) == groups.to.compare)
    selection <- WC84_Fst[row.choice,] > Fst.cutoff & !is.na(WC84_Fst[row.choice, ])
    SNP_genotypes <- combo.NApass[,c(rep(TRUE, times=num_loc_cols), selection)]
    title.text <- paste0("PCA of ",region.text,": genotypes Fst>", Fst.cutoff," loci between ", groups.to.compare, sep="")
  } 	else {
    mean_genotypes <- colMeans(combo.NApass[(num_loc_cols+1):length(combo.NApass[1,])], na.rm=TRUE)
    SNP_genotypes <- combo.NApass[,c(rep(TRUE, times=num_loc_cols), mean_genotypes > 0 & mean_genotypes < 2)]
    title.text <- paste0(str_extract(region.text, "^[:alnum:]*"), sep="")
  }
  # PCA and graph of individuals that passed above filters
  SNP.included <- SNP_genotypes[SNP_genotypes$Fst_group %in% groups.to.plot.PCA,]
  data_for_PCA <- SNP.included[,(num_loc_cols+1):length(SNP.included[1,])]
  impute.pca <- pca(data_for_PCA,method="svdImpute",center=TRUE,scale="none",nPcs=axes)  
  # nipals and svdImpute give similar patterns above
  var_explained <- R2cum(impute.pca)
  scores <- as.data.frame(scores(impute.pca))
  if (flip1 == TRUE){
    scores$PC1 <- -1 * scores$PC1
  }
  if (flip2 == TRUE){
    scores$PC2 <- -1 * scores$PC2
  }
  quartz(title=title.text, 8, 7)
  plot(scores$PC1,scores$PC2, asp=1, pch=23, cex=2, xlab=paste("PC1 (", round(impute.pca@R2[1] * 100, digits=1), "%)", sep=""), ylab=paste("PC2 (", round(impute.pca@R2[2] * 100, digits=1), "%)", sep=""), main=title.text) #EM added cex=2, changed axis labels
  for (i in 1:length(groups.to.plot.PCA)) {
    selection <- SNP.included$Fst_group == groups.to.plot.PCA[i]
    points(scores$PC1[selection], scores$PC2[selection], pch=23, cex=2, bg = group.colors.PCA[i]) #EM added cex
  }
  # return the variance explained, PCA coordinates, and data matrix for all individuals included
  return(list(var_explained=var_explained, scores=scores, data=SNP.included))
}


# Genotype by individual plot
# Graph genotypes of high-Fst loci along a chromosome
# choose only loci that are variable in the dataset (SNPs), and above an Fst threshhold
plotGenotypeByIndividual <- function(groups.to.compare, WC84_Fst.cutoff, missing.fraction.max,
                                     start.pos, end.pos, pos, WC84_Fst, combo.NApass, 
                                     num_loc_cols, freqs, plot.groups, plot.group.colors,
                                     option, group1, group2, chr) {
  
  row.choice <- which(rownames(WC84_Fst) == groups.to.compare)
  selection <- (WC84_Fst[row.choice,] > WC84_Fst.cutoff) & (!is.na(WC84_Fst[row.choice,])) & (pos$position > start.pos) & (pos$position < end.pos)
  SNP.genotypes <- combo.NApass[,c(rep(TRUE, times=num_loc_cols), selection)]
  SNP.freqs <- freqs[, selection]	
  SNP.positions <- pos$position[selection]
  SNP.genotypes.subset <- SNP.genotypes[SNP.genotypes$Fst_group %in% plot.groups,]
  if (option == 2) {
    alt.allele.hi.in.group1 <- which(SNP.freqs[rownames(SNP.freqs)==group1,] > 0.5) + num_loc_cols
    SNP.genotypes.subset[,alt.allele.hi.in.group1] <- -1*SNP.genotypes.subset[,alt.allele.hi.in.group1] + 2
  }
  # Choose sorting order (by group defined above, or by plot_order column in input metadata file)
  #sorted.SNP.genotypes.subset <- SNP.genotypes.subset[order(SNP.genotypes.subset$group, SNP.genotypes.subset$ID),]
  sorted.SNP.genotypes.subset <- SNP.genotypes.subset[order(SNP.genotypes.subset$plot_order, decreasing=TRUE),]
  num.inds <- length(sorted.SNP.genotypes.subset[,1])   # numbering the individuals
  
  # filter out the SNPs that have too much missing data:
  missing.fraction <- colSums(is.na(sorted.SNP.genotypes.subset)) / num.inds
  selection <- missing.fraction <= missing.fraction.max
  sorted.SNP.genotypes.subset.NApass <- sorted.SNP.genotypes.subset[, selection]
  selection2 <- selection[(num_loc_cols+1):length(selection)]
  SNP.positions.NApass <- SNP.positions[selection2]
  
  # Set up the plot window:
  quartz(title=paste0(region.text,": genotypes Fst>", WC84_Fst.cutoff," loci between ", groups.to.compare, sep=""), width=10, height=12)
  chr.length <- max(pos$position)
  genotype.colors <- c("#3f007d", "#807dba", "#dadaeb", "grey50")  # purple shades from colorbrewer
  # choose which type to plot: 1=nucleotide positions; 2=spaced evenly; 3= both
  plot.along.chromosome.type <- 2   
  # plot along chromosome by nucleotide position:
  if (plot.along.chromosome.type == 1 | plot.along.chromosome.type == 3) {
    plot(x=NULL, y=NULL, xlim=c(start.pos, end.pos+0.05*(end.pos-start.pos)), ylim=c(0, num.inds+1), main=paste0("Chr ",chr,": genotypes for Fst>", WC84_Fst.cutoff," loci between ", groups.to.compare, sep=""),xlab=paste0("Location along chromosome ",chr, sep=""), ylab="Individual")
    #plot(x=NULL, y=NULL, xlim=c(850000, 910000), ylim=c(0, num.inds+1))
    # generate my own plotting symbol (a rectangle)
    symbol.x <- c(-0.1, -0.1, 0.1, 0.1, -0.1)
    symbol.y <- c(1, -1, -1, 1, 1)
    plot.symbol <- cbind(symbol.x, symbol.y)
    symbol.size.x <- 5000  #width of box in nucleotides
    symbol.size.y <- 0.8  #height of box in units of individuals
    # cycle through individuals, graphing each type of genotype:
    for (i in 1:num.inds) {
      y <- (-1*i)+num.inds+1  #reverses order of plot top-bottom
      lines(x = c(start.pos,end.pos), y = c(y,y), col = "grey")
      text(x = end.pos, y = y, labels=unlist(strsplit(as.character(sorted.SNP.genotypes.subset.NApass$ID[i]), split='_', fixed=TRUE))[3], cex=0.3, pos=4)
      genotypes <- sorted.SNP.genotypes.subset.NApass[i, (num_loc_cols+1):length(sorted.SNP.genotypes.subset.NApass[1,])]
      hom.ref.locs <- SNP.positions.NApass[genotypes == 0 & !is.na(genotypes)]
      my.symbols(hom.ref.locs, rep(y, times=length(hom.ref.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="red")
      het.locs <- SNP.positions.NApass[genotypes == 1 & !is.na(genotypes)]
      my.symbols(het.locs, rep(y, times=length(het.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="orange")
      hom.alt.locs <- SNP.positions.NApass[genotypes == 2 & !is.na(genotypes)]
      my.symbols(hom.alt.locs, rep(y, times=length(hom.alt.locs)), plot.symbol, xsize=symbol.size.x, ysize=symbol.size.y, col="yellow")
    }
  }
  # plot evenly spaced by SNP order along chromosome:
  # make top part of fig (genotypes for individuals)
  if (plot.along.chromosome.type == 2 | plot.along.chromosome.type == 3) {
    num.SNPs.to.plot <- length(SNP.positions.NApass)
    plot(x=NULL, y=NULL, xlim=c(0.5-0.07*(num.SNPs.to.plot+0.5), 1.07*(num.SNPs.to.plot+0.5)), ylim=c(0.5-0.25*num.inds, num.inds+1), main=paste0(region.text,": genotypes for WC84_Fst>", WC84_Fst.cutoff," loci between ", groups.to.compare, sep=""),xlab=paste0("Order along chromosome ",chr, sep=""), ylab="Individual")
    image.matrix <- t(as.matrix(sorted.SNP.genotypes.subset.NApass[, (num_loc_cols+1):length(sorted.SNP.genotypes.subset.NApass[1,])]))
    image.matrix[is.na(image.matrix)] <- 3

    group.color.box.loc.right <- 1.055*(num.SNPs.to.plot+0.5)
    group.color.box.loc.left <- 0.5-0.055*(num.SNPs.to.plot+0.5)
    box.width <- 0.005*num.SNPs.to.plot * 2
    group.color.box.x.right <- c(-box.width, -box.width, box.width, box.width, -box.width) + group.color.box.loc.right
    group.color.box.x.left <- c(-box.width, -box.width, box.width, box.width, -box.width) + group.color.box.loc.left
    group.color.box.y <- c(0.4, -0.4, -0.4, 0.4, 0.4)
    
    for (i in 1:num.inds) {
      y <- (-1*i)+num.inds+1  #reverses order of plot top-bottom
      name.text.bits <- unlist(strsplit(as.character(sorted.SNP.genotypes.subset.NApass$ID[y]), split='_', fixed=TRUE))
      label.text <- name.text.bits[length(name.text.bits)]
      text(x = num.SNPs.to.plot+0.5, y = y, labels=label.text, cex=0.3, pos=4)
      text(x = 0.5, y = y, labels=label.text, cex=0.3, pos=2)
      polygon(group.color.box.x.right, y+group.color.box.y, border=NA, col=plot.group.colors[which(plot.groups==sorted.SNP.genotypes.subset.NApass$Fst_group[y])])
      polygon(group.color.box.x.left, y+group.color.box.y, border=NA, col=plot.group.colors[which(plot.groups==sorted.SNP.genotypes.subset.NApass$Fst_group[y])])
    }
    
    # generate my own plotting symbol (a rectangle)
    symbol.x <- c(-0.5, -0.5, 0.5, 0.5, -0.5)
    symbol.y <- c(0.4, -0.4, -0.4, 0.4, 0.4)
    # generate triangles for plotting heterozygotes
    triangle1.x <- c(-0.5, -0.5, 0.5, -0.5)
    triangle1.y <- c(0.4, -0.4, 0.4, 0.4)
    triangle2.x <- c(-0.5, 0.5, 0.5, -0.5)
    triangle2.y <- c(-0.4, -0.4, 0.4, -0.4)
    # cycle through individuals, graphing each type of genotype:
    for (i in 1:num.inds) {
      y <- (-1*i)+num.inds+1  #reverses order of plot top-bottom
      lines(x = c(0.5,num.SNPs.to.plot+0.5), y = c(y,y), col = "grey40")
      genotypes <- sorted.SNP.genotypes.subset.NApass[y, (num_loc_cols+1):length(sorted.SNP.genotypes.subset.NApass[1,])]
      hom.ref.locs <- which(genotypes == 0 & !is.na(genotypes))
      if (length(hom.ref.locs) > 0) {
        for (j in 1:length(hom.ref.locs)) {
          polygon(hom.ref.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[1])
        }
      }
      het.locs <- which(genotypes == 1 & !is.na(genotypes))
      if (length(het.locs) > 0) {
        for (j in 1:length(het.locs)) {
          #polygon(het.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[2])  # draws rectangle in hetero color
          polygon(het.locs[j]+triangle1.x, y+triangle1.y, border=NA, col=genotype.colors[1])  # draws triangle in hom ref color
          polygon(het.locs[j]+triangle2.x, y+triangle2.y, border=NA, col=genotype.colors[3])  # draws triangle in hom alt color
        }
      }
      hom.alt.locs <- which(genotypes == 2 & !is.na(genotypes))
      if (length(hom.alt.locs) > 0) {
        for (j in 1:length(hom.alt.locs)) {
          polygon(hom.alt.locs[j]+symbol.x, y+symbol.y, border=NA, col=genotype.colors[3])
        }
      }
    }
  }
  # make lower part of figure (indicating position along chromosome) #NOTE THAT CHROMOSOME LENGTH NOT QUITE THE TRUE LENGTH
  chr.line.y <- 0.5-0.2*num.inds
  top.hatch.line.y1 <- 0.5-0.005*num.inds
  top.hatch.line.y2 <- 0.5-0.02*num.inds
  low.hatch.line.y1 <- 0.5-0.18*num.inds
  low.hatch.line.y2 <- 0.5-0.2*num.inds
  lines(x = c(0.5,num.SNPs.to.plot+0.5), y = c(chr.line.y,chr.line.y), lwd=4, col = "black") #draws chromosome line
  text(x=(0.5+(num.SNPs.to.plot+0.5)/2), y=chr.line.y-0.05*num.inds, paste0("Location along chromosome ",chr, sep=""))
  text(x=0.5, y=chr.line.y-0.025*num.inds, start.pos)
  text(x=num.SNPs.to.plot+0.5, y=chr.line.y-0.025*num.inds, end.pos)
  chr.plot.ratio <- num.SNPs.to.plot/(end.pos-start.pos)
  for (i in 1:length(SNP.positions.NApass)) {
    lines(x=c(i,i), y=c(top.hatch.line.y1, top.hatch.line.y2), lwd=0.5, col="grey20")
    lines(x=c(i, 1+chr.plot.ratio*(SNP.positions.NApass[i]-start.pos)), y=c(top.hatch.line.y2, low.hatch.line.y1), lwd=0.5, col="grey20")
    lines(x=c(1+chr.plot.ratio*(SNP.positions.NApass[i]-start.pos), 1+chr.plot.ratio*(SNP.positions.NApass[i]-start.pos)), y=c(low.hatch.line.y1, low.hatch.line.y2), lwd=0.5, col="grey20")
  }
}


# plot Genome-wide plot of windowed Fst, Dxy, pi
# read files for each chromosome, and plot Fst, Dxy, and pi for each chromosome in a single window
plotGenomeFstDxyPi <- function(base.name, tag.name, window_size, chromosomes.to.plot, 
                               max_Dxy_axis, max_pi_axis, transparency, line_transparency,
                               groups.to.plot.WC84_Fst, group.colors.WC84_Fst,
                               groups.to.plot.Dxy, group.colors.Dxy,
                               groups.to.plot.pi, group.colors.pi) {
  #plot all chromosomes in one quartz window
  # plot dimensions in inches:
  window.width <- 12
  window.height <- 12
  # numbers below are in inches:
  left.margin <- 1 / window.width  # the left outer margin of the figure (number is inches)
  right.margin <- 1 / window.width
  top.margin <- 0.5 / window.height
  bottom.margin <- 0.5 / window.height
  row.height <- 1 / window.height  # height of each row of the figure (row contains Fst, Dxy, and pi plot)
  plot.height <- row.height / 3
  row.space <- 0.25 / window.height  # gap between rows 
  bp.per.inch <- 20000000  # number of bp of sequence per inch
  gap.between.chr <- 0.6 # gap in inches between chromosomes plots
  gap.between.chr.bp <- gap.between.chr*bp.per.inch
  
  # get sizes of chromosomes (note this actually isn't the true length, just the mean position of the rightmost window):
  chr.length <- NULL	
  for (i in 1:length(chromosomes.to.plot)) {
    chr.text <- paste0("Chr",chromosomes.to.plot[i], "_whole", sep="")
    #load(paste0("~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_SNP_summary_stats_by_chromosome_from_R/GW_SNP_RollingMean_stats_GW_Lane5_plus_Liz.45samples.troch_vir_plumb.",load.text,"_from_R.R"))
    load(paste0(base.name,tag.name,chr.text,"_window",window_size,"_WindowStats_from_R.R"))
    #chr.length[i] <- rolling.mean.pos.Fst[length(rolling.mean.pos.Fst)]
    chr.length[i] <- rolling.mean.pos.WC84_Fst[length(rolling.mean.pos.WC84_Fst)]
  }
  
  chr.plotted.already <- rep(FALSE, length(chromosomes.to.plot))
  
  # make plot order by row
  inch.per.row <- (window.width-(left.margin+right.margin)*window.width)
  bp.per.row <- bp.per.inch * inch.per.row
  row <- 1
  chr.in.row <- matrix(NA, 20, 20)
  chr.bp.start.in.row <- matrix(NA, 20, 20)
  bp.start <- 0
  remaining.row.length <- bp.per.row
  while (sum(chr.plotted.already==FALSE) > 0) {  # repeat until all chr have a row to be plotted in
    place.in.row <- 1
    # look through chromosomes
    for (i in 1:length(chromosomes.to.plot)) {
      # if not plotted and short enough, add to row:
      if (chr.plotted.already[i]==FALSE && chr.length[i]<=remaining.row.length) {
        chr.in.row[row, place.in.row] <- chromosomes.to.plot[i]
        chr.bp.start.in.row[row, place.in.row] <- bp.start
        remaining.row.length <- remaining.row.length - chr.length[i] - gap.between.chr.bp
        chr.plotted.already[i] <- TRUE
        place.in.row <- place.in.row+1
        bp.start <- bp.start + chr.length[i] + gap.between.chr.bp
      }		
    }
    row <- row+1
    bp.start <- 0
    remaining.row.length <- bp.per.row
  }
  
  # Plot in the order defined above
  
  quartz(title="WC84_Fst, Dxy, and pi for all chromosomes", width=window.width, height=window.height)  # width and height in inches
  # plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt='n', yaxt='n', xlab=NULL, ylab= NULL, bty='n', las=1)
  
  num.rows <- sum(!is.na(chr.in.row[,1]))
  new.fig.parameter <- FALSE
  for (row.num in 1:num.rows) {
    num.chr.in.row <- sum(!is.na(chr.in.row[row.num,]))
    for (order.num in 1:num.chr.in.row) {
      chr.text <- paste0("Chr",chr.in.row[row.num,order.num], "_whole", sep="")
      #load(paste0("~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_SNP_summary_stats_by_chromosome_from_R/GW_SNP_RollingMean_stats_",load.text,"_from_R.R"))
      load(paste0(base.name,tag.name,chr.text,"_window",window_size,"_WindowStats_from_R.R"))
      print(paste0("Loaded saved rolling mean stats for Chr ", chr.in.row[row.num,order.num], sep=""))
      plot.bp.length <- chr.length[chromosomes.to.plot==chr.in.row[row.num,order.num]]
      # plot WC84_Fst:
      fig.left.loc <- left.margin + (chr.bp.start.in.row[row.num,order.num]/bp.per.row)*(inch.per.row/window.width)   
      fig.right.loc <- fig.left.loc + (plot.bp.length/bp.per.row)*(inch.per.row/window.width)
      par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+plot.height+(row.num-1)*(row.height+row.space)), 1-(top.margin+(row.num-1)*(row.height+row.space))), new = new.fig.parameter, mai=c(0,0,0,0), cex=0.5)  # mai sets the margins of the plot to zero, matching
      new.fig.parameter <- TRUE
      #upper.xlim <- bp.per.inch*(window.width-(left.margin+right.margin)*window.width)
      plot(x=NULL, y=NULL, xlim=c(0, plot.bp.length), ylim=c(0, 1.2), yaxp=c(0, 1, n=1), ylab=expression(italic('F')[ST]*"  "), xlab=NA, xaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex=1, cex.lab=1.55, cex.axis=0.75)
      text(-2000000, 1.4, chr.in.row[row.num,order.num], pos=4, cex=2, xpd=NA)
      for (j in 1:length(groups.to.plot.WC84_Fst)) {
        color.rgb <- col2rgb(group.colors.WC84_Fst[j]) /255   # divide to convert color scales to 0-1
        lines(rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst[rownames(rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[j],],
              col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line_transparency))
        xx <- c(rolling.mean.pos.WC84_Fst[1], rolling.mean.pos.WC84_Fst, rolling.mean.pos.WC84_Fst[length(rolling.mean.pos.WC84_Fst)])
        yy <- c(0, rolling.mean.WC84_Fst[rownames(rolling.mean.WC84_Fst)==groups.to.plot.WC84_Fst[j],], 0)
        polygon(xx, yy, col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=transparency), border=NA)
      }
      # plot Dxy:
      par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+2*plot.height+(row.num-1)*(row.height+row.space)), 1-(top.margin+plot.height+(row.num-1)*(row.height+row.space))), new = TRUE, mai=c(0,0,0,0))  # mai sets the margins of the plot to zero, matching
      plot(x=NULL, y=NULL, xlim=c(0, plot.bp.length), ylim=c(0, max_Dxy_axis), yaxp=c(0, max_Dxy_axis*5/8, n=1), ylab=expression(italic('D')[xy]*"    "), xlab=NA, xaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex.lab=1.5, cex.axis=0.75)  # 
      for (j in 1:length(groups.to.plot.Dxy)) {
        color.rgb <- col2rgb(group.colors.Dxy[j]) /255   # divide to convert color scales to 0-1
        lines(rolling.mean.pos.Dxy, rolling.mean.Dxy[rownames(rolling.mean.Dxy)==groups.to.plot.Dxy[j],], 
              col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line_transparency))
        xx <- c(rolling.mean.pos.Dxy[1], rolling.mean.pos.Dxy, rolling.mean.pos.Dxy[length(rolling.mean.pos.Dxy)])
        yy <- c(0, rolling.mean.Dxy[rownames(rolling.mean.Dxy)==groups.to.plot.Dxy[j],], 0)
        polygon(xx, yy, col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=transparency), border=NA)
      }
      # plot pi: 
      par(fig=c(fig.left.loc, fig.right.loc, 1-(top.margin+3*plot.height+(row.num-1)*(row.height+row.space)), 1-(top.margin+2*plot.height+(row.num-1)*(row.height+row.space))), new = TRUE, mai=c(0,0,0,0))  # mai sets the margins of the plot to zero, matching
      plot(x=NULL, y=NULL, xlim=c(0, plot.bp.length), ylim=c(0, max_pi_axis), yaxp=c(0, max_pi_axis*5/8, n=1), ylab=expression(pi*"    "), xlab=NA, xaxt='n', bty='n', las=1, xaxs="i", tcl=-0.25, xpd=NA, cex.lab=1.5, cex.axis=0.75)
      title(xlab="Location", outer=TRUE)
      for (j in 1:length(groups.to.plot.pi)) {
        color.rgb <- col2rgb(group.colors.pi[j]) /255   # divide to convert color scales to 0-1
        lines(rolling.mean.pos.pi, rolling.mean.pi[rownames(rolling.mean.pi)==groups.to.plot.pi[j],], 
              col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=line_transparency))
        xx <- c(rolling.mean.pos.pi[1], rolling.mean.pos.pi, rolling.mean.pos.pi[length(rolling.mean.pos.pi)])
        yy <- c(0, rolling.mean.pi[rownames(rolling.mean.pi)==groups.to.plot.pi[j],], 0)
        polygon(xx, yy, col=rgb(color.rgb[1], color.rgb[2], color.rgb[3], alpha=transparency), border=NA)
      }
    }
  }
}                

# compile windowed WC84_Fst, Dxy, pi for a bunch of chromosomes,
# using previously saved files
compileWindowedStats <- function(base.name, tag.name, chromosomes.to.combine, window.size) {
  genome.rolling.mean.positions <- NULL
  genome.rolling.mean.WC84_Fst <- NULL
  genome.rolling.mean.Dxy <- NULL
  genome.rolling.mean.pi <- NULL
  for (i in 1:length(chromosomes.to.combine)) {
    chr.text <- paste0("Chr",chromosomes.to.combine[i], "_whole", sep="")
    load(paste0(base.name,tag.name,chr.text,"_window",window_size,"_WindowStats_from_R.R"))
    chromosome.ID.vector <- rep(chromosomes.to.combine[i], times=length(rolling.mean.pos.WC84_Fst))
    genome.rolling.mean.positions <- cbind(genome.rolling.mean.positions, rbind(chromosome.ID.vector, rolling.mean.pos.WC84_Fst))
    genome.rolling.mean.WC84_Fst <- cbind(genome.rolling.mean.WC84_Fst, rolling.mean.WC84_Fst)
    genome.rolling.mean.Dxy <- cbind(genome.rolling.mean.Dxy, rolling.mean.Dxy)
    genome.rolling.mean.pi <- cbind(genome.rolling.mean.pi, rolling.mean.pi)
  }
  return(list(positions=genome.rolling.mean.positions, 
              WC84_Fst=genome.rolling.mean.WC84_Fst, 
              Dxy=genome.rolling.mean.Dxy, 
              pi=genome.rolling.mean.pi))
}


# For compiled windowed data ----

# Make a plot of windowed Fst vs. Dxy 
# using compiled genome-wide info, and conduct pearson's correlation test:
plotFst_Dxy <- function(comp, autosome.genome.rolling.stats, method) {
  quartz(title=paste0("WC84_Fst vs. Dxy, based on all autosomes, ", comp, sep=""), width=4, height=4)
  row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comp)
  plot(autosome.genome.rolling.stats$WC84_Fst[row.choice,], autosome.genome.rolling.stats$Dxy[row.choice,], 
       pch=16, cex=0.25, xlab=paste0("windowed WC84_Fst for ", comp, sep=""),
       ylab=paste0("windowed Dxy for ", comp, sep=""))
  result <- cor.test(autosome.genome.rolling.stats$WC84_Fst[row.choice,], autosome.genome.rolling.stats$Dxy[row.choice,],
                     method=method)
  return(result)
}


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


# make scatterplots of Fst vs. Dxy for three group comparisons.
# each of the "comp_" objects should be have two group names and a color name, e.g.:
# comp1 <- c("troch", "vir", "green3") 
# returns the statistical correlation tests, using the cor.method as specified
plots3Fst_Dxy <- function(comp1, comp2, comp3,
                          autosome.genome.rolling.stats, cor.method) {
  quartz(title=paste0("Scatterplots of Fst vs. Dxy in three comparisons"), width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.5
  Fst.row.choice1 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp1[1], comp1[2], sep="_"))
  Fst.row.choice2 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp2[1], comp2[2], sep="_"))
  Fst.row.choice3 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp3[1], comp3[2], sep="_"))
  Fst.vector1 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice1,]
  Fst.vector2 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice2,]
  Fst.vector3 <- autosome.genome.rolling.stats$WC84_Fst[Fst.row.choice3,]
  Dxy.row.choice1 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp1[1], comp1[2], sep="_"))
  Dxy.row.choice2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp2[1], comp2[2], sep="_"))
  Dxy.row.choice3 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp3[1], comp3[2], sep="_"))
  Dxy.vector1 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice1,]
  Dxy.vector2 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice2,]
  Dxy.vector3 <- autosome.genome.rolling.stats$Dxy[Dxy.row.choice3,]
  high.Dxy <- round(max(c(Dxy.vector1, Dxy.vector2, Dxy.vector3)), 3) + 0.001
  # first plot:
  plot(x=Fst.vector1, y=Dxy.vector1, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  title(ylab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  # Add cubic spline:
  lines(smooth.spline(x=Fst.vector1, y=Dxy.vector1, spar=1), lty = 1, col = alpha(comp1[3], 0.75), lwd=2)  # 0.8 is transparency
  # second plot:
  plot(x=Fst.vector2, y=Dxy.vector2, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  title(ylab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  lines(smooth.spline(x=Fst.vector2, y=Dxy.vector2, spar=1), lty = 1, col = alpha(comp2[3], 0.75), lwd=2)
  # third plot:
  plot(x=Fst.vector3, y=Dxy.vector3, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  title(ylab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  lines(smooth.spline(x=Fst.vector3, y=Dxy.vector3, spar=1), lty = 1, col = alpha(comp3[3], 0.75), lwd=2)
  # do statistical tests of correlation:
  test1 <- cor.test(Fst.vector1, Dxy.vector1, method=cor.method)
  test2 <- cor.test(Fst.vector2, Dxy.vector2, method=cor.method)
  test3 <- cor.test(Fst.vector3, Dxy.vector3, method=cor.method)
  return(list(test1=test1, test2=test2, test3=test3))
}


# Fst_Fst plot
# Make a bivariate plot of WC84_Fst using for two population pair comparisons:
# return the correlation, according to cor.method
plotFst_Fst <- function(comp1, comp2, cor.method, 
                        autosome.genome.rolling.stats) {
  quartz(title=paste0("Bivariate plot of windowed WC84_Fst, based on all autosomes, ", comp1," vs ", comp2, sep=""), width=6, height=6)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comp1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comp2)
  plot(autosome.genome.rolling.stats$WC84_Fst[row.choice.1,], autosome.genome.rolling.stats$WC84_Fst[row.choice.2,], cex=0.5)
  # lines(c(0,1), c(0,1))
  test <- cor.test(autosome.genome.rolling.stats$WC84_Fst[row.choice.1,], autosome.genome.rolling.stats$WC84_Fst[row.choice.2,], method=cor.method)
  return(test)
}


# Make a bivariate plot of pi in one group and Dxy for one pair,
# and return statistical test of correlation, using method cor.method:
plotPi_Dxy <- function(group1, comp2, cor.method, 
                       autosome.genome.rolling.stats) {
  quartz(title=paste0("Bivariate plot of pi in ", group1," vs. Dxy between", comp2, sep=""), width=6, height=6)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == comp2)
  plot(autosome.genome.rolling.stats$pi[row.choice.1,], autosome.genome.rolling.stats$Dxy[row.choice.2,], cex=0.5)
  test <- cor.test(autosome.genome.rolling.stats$pi[row.choice.1,], autosome.genome.rolling.stats$Dxy[row.choice.2,])
  return(test)
}


# Fst_Fst 3 scatterplot matrix
# make scatterplots of Fst vs. Fst for three comparisons of paired Fst
# returns the statistical correlation tests, using the cor.method as specified
FstMatrixPlot <- function(comp1, comp2, comp3, cor.method, 
                          autosome.genome.rolling.stats) {
  quartz(title=paste0("Scatterplot matrix of windowed WC84_Fst, based on all autosomes, ", comp1," vs ", comp2, " vs ", comp3, sep=""), width=6, height=6)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comp1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comp2)
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == comp3)
  vector1 <- as.vector(autosome.genome.rolling.stats$WC84_Fst[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$WC84_Fst[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$WC84_Fst[row.choice.3,])
  pairs(cbind(vector1, vector2, vector3), labels=c(comp1, comp2, comp3), pch=16, cex=0.25,
        lower.panel=NULL, diag.panel=NULL)
  test1 <- cor.test(vector1, vector2, method=cor.method)
  test2 <- cor.test(vector1, vector3, method=cor.method)
  test3 <- cor.test(vector2, vector3, method=cor.method)
  return(list(test1=test1, test2=test2, test3=test3))
}


# make scatterplots of Fst vs. Fst for three comparisons of group pairs
# each of the "comp_" objects should be have two group names and a color name, e.g.:
# comp1 <- c("troch", "vir", "green3") 
# returns the statistical correlation tests, using the cor.method as specified
plots3Fst_Fst <- function(comp1, comp2, comp3, 
                          autosome.genome.rolling.stats, cor.method) {
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp1[1], comp1[2], sep="_"))
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp2[1], comp2[2], sep="_"))
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == paste(comp3[1], comp3[2], sep="_"))
  vector1 <- as.vector(autosome.genome.rolling.stats$WC84_Fst[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$WC84_Fst[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$WC84_Fst[row.choice.3,])
  quartz(title=paste0("Scatterplots of Fst in three comparisons"), width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.5
  # first plot:
  plot(x=vector1, y=vector3, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,1), yaxp=c(0,1,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp1[2])) * " to " * italic(.(comp1[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp1[2])) * " to " * italic(.(comp1[1]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  title(ylab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2, cex.lab=label.size, col.lab=comp3[3])
  # second plot:
  plot(x=vector2, y=vector1, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,1), yaxp=c(0,1,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  title(ylab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2, cex.lab=label.size, col.lab=comp1[3])
  # third plot:
  plot(x=vector3, y=vector2, xlim=c(0,1), xaxp=c(0,1,2), ylim=c(0,1), yaxp=c(0,1,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp3[2])) * " to " * italic(.(comp3[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp3[2])) * " to " * italic(.(comp3[1]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  title(ylab=bquote(italic('F')[ST] * ", " * phantom(italic(.(comp2[2])) * " to " * italic(.(comp2[1])))), line=2, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('F')[ST] * ", ") * italic(.(comp2[2])) * " to " * italic(.(comp2[1]))), line=2, cex.lab=label.size, col.lab=comp2[3])
  # statistical tests
  test1 <- cor.test(vector1, vector3, method=cor.method)
  test2 <- cor.test(vector2, vector1, method=cor.method)
  test3 <- cor.test(vector3, vector2, method=cor.method)
  return(list(test1=test1, test2=test2, test3=test3))
}

# make scatterplots of Dxy vs. Dxy for three comparisons of group pairs
# each of the "comp_" objects should be have two group names and a color name, e.g.:
# comp1 <- c("troch", "vir", "green3") 
# returns the statistical correlation tests, using the cor.method as specified
plots3Dxy_Dxy <- function(comp1, comp2, comp3, 
                          autosome.genome.rolling.stats, cor.method) {
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp1[1], comp1[2], sep="_"))
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp2[1], comp2[2], sep="_"))
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$Dxy) == paste(comp3[1], comp3[2], sep="_"))
  vector1 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.3,])
  high.Dxy <- round(max(c(vector1, vector2, vector3)), 3) + 0.00
  quartz(title=paste0("Scatterplots of Dxy in three comparisons"), width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.5
  # first plot:
  plot(x=vector1, y=vector3, xlim=c(0,high.Dxy), xaxp=c(0,high.Dxy,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp1[2])) * " to " * italic(.(comp1[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp1[2])) * " to " * italic(.(comp1[1]))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  title(ylab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp3[1])) * " to " * italic(.(comp3[2])))), line=2, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp3[1])) * " to " * italic(.(comp3[2]))), line=2, cex.lab=label.size, col.lab=comp3[3])
  # second plot:
  plot(x=vector2, y=vector1, xlim=c(0,high.Dxy), xaxp=c(0,high.Dxy,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp2[1])) * " to " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp2[1])) * " to " * italic(.(comp2[2]))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  title(ylab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp1[1])) * " to " * italic(.(comp1[2])))), line=2, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp1[1])) * " to " * italic(.(comp1[2]))), line=2, cex.lab=label.size, col.lab=comp1[3])
  # third plot:
  plot(x=vector3, y=vector2, xlim=c(0,high.Dxy), xaxp=c(0,high.Dxy,2), ylim=c(0,high.Dxy), yaxp=c(0,high.Dxy,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp3[2])) * " to " * italic(.(comp3[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp3[2])) * " to " * italic(.(comp3[1]))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  title(ylab=bquote(italic('D')[xy] * ", " * phantom(italic(.(comp2[2])) * " to " * italic(.(comp2[1])))), line=2, cex.lab=label.size)
  title(ylab=bquote(phantom(italic('D')[xy] * ", ") * italic(.(comp2[2])) * " to " * italic(.(comp2[1]))), line=2, cex.lab=label.size, col.lab=comp2[3])
  # statistical tests
  test1 <- cor.test(vector1, vector3, method=cor.method)
  test2 <- cor.test(vector2, vector1, method=cor.method)
  test3 <- cor.test(vector3, vector2, method=cor.method)
  return(list(test1=test1, test2=test2, test3=test3))
}



# Make a scatterplot matrix for Dxy using above compiled genome-wide info.
DxyMatrixPlot <- function(comp1, comp2, comp3, cor.method, 
                          autosome.genome.rolling.stats, ...) {
  quartz(title=paste0("Scatterplot matrix of windowed Dxy, based on compiled genome, ", comp1," vs ", comp2, " vs ", comp3, sep=""), width=6, height=6)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$Dxy) == comp1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == comp2)
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$Dxy) == comp3)
  vector1 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.3,])
  pairs(cbind(vector1, vector2, vector3), labels=c(comp1, comp2, comp3), pch=23, cex=0.25,
        ...)
  test1 <- cor.test(vector1, vector2, method=cor.method)
  test2 <- cor.test(vector1, vector3, method=cor.method)
  test3 <- cor.test(vector2, vector3, method=cor.method)
  return(list(test1=test1, test2=test2, test3=test3))
}


# RND hist
# Make a histogram of relative node depth 
# (Dxy for one pair of pops standardized by Dxy to an outgroup).
# Example comparisons:
# comp1 <- "troch_plumb"   # the focal comparison
# comp2 <- "troch_vir"    # taxon A to the outgroup
# comp3 <- "vir_plumb"  # taxon B to the outgroup
RNDhist <- function(comp1, comp2, comp3, autosome.genome.rolling.stats) {
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$Dxy) == comp1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == comp2)
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$Dxy) == comp3)
  vector1 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$Dxy[row.choice.3,])
  ave.Dxy.to.outgroup <- vector2 + vector3 / 2
  RND <- vector1 / ave.Dxy.to.outgroup
  quartz(title=paste0("Histogram of RND: ", comp1," standardized by ", comp2, " and ", comp3, sep=""), width=6, height=6)
  hist(RND, 30)
  return(RND)
}


# Pi scatterplot matrix 3groups
PiMatrixPlot3 <- function(group1, group2, group3, cor.method, 
                          autosome.genome.rolling.stats) {
  quartz(title=paste0("Scatterplot matrix of windowed pi, based on compiled chromosomes, ", group1," vs ", group2, " vs ", group3, sep=""), width=6, height=6)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$pi) == group3)
  vector1 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.3,])
  pairs(cbind(vector1, vector2, vector3), labels=c(group1, group2, group3), pch=16, cex=0.25,
        lower.panel=NULL, diag.panel=NULL)
  test1 <- cor.test(vector1, vector2, method=cor.method)
  test2 <- cor.test(vector1, vector3, method=cor.method)
  test3 <- cor.test(vector2, vector3, method=cor.method) 
  return(list(test1=test1, test2=test2, test3=test3))
}


# Pi scatterplot matrix 4groups
PiMatrixPlot4 <- function(group1, group2, group3, group4, cor.method, 
                          autosome.genome.rolling.stats, ...) {
  quartz(title=paste0("Scatterplot matrix of windowed pi, based on compiled chromosomes, ", group1," vs ", group2, " vs ", group3, " vs ", group4, sep=""), width=6, height=6)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$pi) == group3)
  row.choice.4 <- which(rownames(autosome.genome.rolling.stats$pi) == group4)
  vector1 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.3,])
  vector4 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.4,])
  pairs(cbind(vector1, vector2, vector3, vector4), labels=c(group1, group2, group3, group4), pch=16, cex=0.25,
        lower.panel=NULL, diag.panel=NULL, ...)
  test1 <- cor.test(vector1, vector2, method=cor.method)
  test2 <- cor.test(vector1, vector3, method=cor.method)
  test3 <- cor.test(vector1, vector4, method=cor.method)
  test4 <- cor.test(vector2, vector3, method=cor.method)
  test5 <- cor.test(vector2, vector4, method=cor.method) 
  test6 <- cor.test(vector3, vector4, method=cor.method) 
  return(list(test1=test1, test2=test2, test3=test3, test4=test4, test5=test5, test6=test6))
}


# Make 3 scatterplots of pi in a line
plots3pi_pi <- function(groups, colors, cor.method, 
                        autosome.genome.rolling.stats) {
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == groups[1])
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == groups[2])
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$pi) == groups[3])
  vector1 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2,])
  vector3 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.3,])
  quartz(title=paste0("Scatterplots of pi in three comparisons"), width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.5
  high.pi <- round(max(c(vector1, vector2, vector3)), 3) + 0.001
  # first plot:
  plot(x=vector1, y=vector2, xlim=c(0,high.pi), xaxp=c(0,high.pi,2), ylim=c(0,high.pi), yaxp=c(0,high.pi,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic(pi) * ", " * phantom(italic(.(groups[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic(pi) * ", ") * italic(.(groups[1]))), line=2.5, cex.lab=label.size, col.lab=colors[1])
  title(ylab=bquote(italic(pi) * ", " * phantom(italic(.(groups[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi) * ", ") * italic(.(groups[2]))), line=2.5, cex.lab=label.size, col.lab=colors[2])
  # second plot:
  plot(x=vector1, y=vector3, xlim=c(0,high.pi), xaxp=c(0,high.pi,2), ylim=c(0,high.pi), yaxp=c(0,high.pi,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic(pi) * ", " * phantom(italic(.(groups[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic(pi) * ", ") * italic(.(groups[1]))), line=2.5, cex.lab=label.size, col.lab=colors[1])
  title(ylab=bquote(italic(pi) * ", " * phantom(italic(.(groups[3])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi) * ", ") * italic(.(groups[3]))), line=2.5, cex.lab=label.size, col.lab=colors[3])
  # third plot:
  plot(x=vector2, y=vector3, xlim=c(0,high.pi), xaxp=c(0,high.pi,2), ylim=c(0,high.pi), yaxp=c(0,high.pi,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic(pi) * ", " * phantom(italic(.(groups[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic(pi) * ", ") * italic(.(groups[2]))), line=2.5, cex.lab=label.size, col.lab=colors[2])
  title(ylab=bquote(italic(pi) * ", " * phantom(italic(.(groups[3])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi) * ", ") * italic(.(groups[3]))), line=2.5, cex.lab=label.size, col.lab=colors[3])
  test1 <- cor.test(vector1, vector2, method=cor.method)
  test2 <- cor.test(vector1, vector3, method=cor.method)
  test3 <- cor.test(vector2, vector3, method=cor.method) 
  return(list(test1=test1, test2=test2, test3=test3))
}


# Make 3 scatterplots of meanPi in a line
plots3meanPi_meanPi <- function(comp1, comp2, comp3, cor.method, 
                        autosome.genome.rolling.stats) {
  row.choice.1A <- which(rownames(autosome.genome.rolling.stats$pi) == comp1[1])
  row.choice.1B <- which(rownames(autosome.genome.rolling.stats$pi) == comp1[2])
  row.choice.2A <- which(rownames(autosome.genome.rolling.stats$pi) == comp2[1])
  row.choice.2B <- which(rownames(autosome.genome.rolling.stats$pi) == comp2[2])
  row.choice.3A <- which(rownames(autosome.genome.rolling.stats$pi) == comp3[1])
  row.choice.3B <- which(rownames(autosome.genome.rolling.stats$pi) == comp3[2])
  vector1A <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1A,])
  vector1B <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1B,])
  vector2A <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2A,])
  vector2B <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2B,])
  vector3A <- as.vector(autosome.genome.rolling.stats$pi[row.choice.3A,])
  vector3B <- as.vector(autosome.genome.rolling.stats$pi[row.choice.3B,])
  vector1 <- (vector1A + vector1B) / 2   # this is mean_Pi
  vector2 <- (vector2A + vector2B) / 2
  vector3 <- (vector3A + vector3B) / 2
  quartz(title=paste0("Scatterplots of meanPi in three comparisons"), width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.3
  high.pi <- round(max(c(vector1, vector2, vector3)), 3) + 0.000
  # first plot:
  plot(x=vector1, y=vector3, xlim=c(0,high.pi), xaxp=c(0,high.pi,2), ylim=c(0,high.pi), yaxp=c(0,high.pi,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote("Mean " * italic(pi) * ", " * phantom(italic(.(comp1[1])) * " & " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom("Mean " * italic(pi) * ", ") * italic(.(comp1[1]) * " & " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  title(ylab=bquote("Mean " * italic(pi) * ", " * phantom(italic(.(comp3[1])) * " & " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom("Mean " * italic(pi) * ", ") * italic(.(comp3[1]) * " & " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  # second plot:
  plot(x=vector2, y=vector1, xlim=c(0,high.pi), xaxp=c(0,high.pi,2), ylim=c(0,high.pi), yaxp=c(0,high.pi,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote("Mean " * italic(pi) * ", " * phantom(italic(.(comp2[1])) * " & " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom("Mean " * italic(pi) * ", ") * italic(.(comp2[1]) * " & " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  title(ylab=bquote("Mean " * italic(pi) * ", " * phantom(italic(.(comp1[1])) * " & " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom("Mean " * italic(pi) * ", ") * italic(.(comp1[1]) * " & " * italic(.(comp1[2])))), line=2.5, cex.lab=label.size, col.lab=comp1[3])
  # third plot:
  plot(x=vector3, y=vector2, xlim=c(0,high.pi), xaxp=c(0,high.pi,2), ylim=c(0,high.pi), yaxp=c(0,high.pi,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote("Mean " * italic(pi) * ", " * phantom(italic(.(comp3[1])) * " & " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom("Mean " * italic(pi) * ", ") * italic(.(comp3[1]) * " & " * italic(.(comp3[2])))), line=2.5, cex.lab=label.size, col.lab=comp3[3])
  title(ylab=bquote("Mean " * italic(pi) * ", " * phantom(italic(.(comp2[1])) * " & " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom("Mean " * italic(pi) * ", ") * italic(.(comp2[1]) * " & " * italic(.(comp2[2])))), line=2.5, cex.lab=label.size, col.lab=comp2[3])
  test1 <- cor.test(vector1, vector3, method=cor.method)
  test2 <- cor.test(vector2, vector1, method=cor.method)
  test3 <- cor.test(vector3, vector2, method=cor.method) 
  return(list(test1=test1, test2=test2, test3=test3))
}

# basic plot of mean_pi vs. Dxy in two taxa:
plotMeanPi_Dxy <- function(group1, group2, cor.method, 
                           autosome.genome.rolling.stats) {
  group.pair <- paste0(group1, "_", group2)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  vector1 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1,])
  vector2 <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2,])
  mean_pi <- (vector1 + vector2) / 2
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair)
  Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  quartz(title=paste0("Scatterplot of mean windowed pi vs. Dxy, based on all autosomes, between", group1,"_", group2, sep=""), width=6, height=6)
  plot(mean_pi, Dxy, cex=0.25)
  test <- cor.test(mean_pi, Dxy, method=cor.method)
  return(test)
}


# Make a scatterplot of mean_pi in two taxa vs. mean_pi in two other taxa.
plotMeanPi_MeanPi <- function(group1A, group1B, group2A, group2B,
                              autosome.genome.rolling.stats, cor.method) {
  row.choice.1A <- which(rownames(autosome.genome.rolling.stats$pi) == group1A)
  row.choice.1B <- which(rownames(autosome.genome.rolling.stats$pi) == group1B)
  row.choice.2A <- which(rownames(autosome.genome.rolling.stats$pi) == group2A)
  row.choice.2B <- which(rownames(autosome.genome.rolling.stats$pi) == group2B)
  vector1A <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1A,])
  vector1B <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1B,])
  mean1 <- (vector1A + vector1B) / 2
  vector2A <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2A,])
  vector2B <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2B,])
  mean2 <- (vector2A + vector2B) / 2
  quartz(title=paste0("Scatterplot of mean windowed pi, based on compiled chromosomes, between ", group1A,"_", group1B, " vs ", group2A,"_", group2B, sep=""), width=6, height=6)
  plot(mean1, mean2, cex=0.25, xlab=NA, ylab=NA)
  test <- cor.test(mean1, mean2)
  title(xlab=bquote("mean " * italic(pi) * " for " * italic(.(group1A)) * " and " * italic(.(group1B))), line=2.5, cex.lab=1)
  title(ylab=bquote("mean " * italic(pi) * " for " * italic(.(group2A)) * " and " * italic(.(group2B))), line=2.5, cex.lab=1)
  return(test)
}


# MeanPi/Dxy corr plot
# Make a scatterplot of mean_pi / Dxy in two taxa vs. mean_pi / Dxy in two other taxa.
plotMeanPiOverDxyCorr <- function(group1A, group1B, color1,
                                  group2A, group2B, color2,
                                  autosome.genome.rolling.stats, cor.method) {
  group.pair_1 <- paste0(group1A, "_", group1B)
  group.pair_2 <- paste0(group2A, "_", group2B)
  row.choice.1A <- which(rownames(autosome.genome.rolling.stats$pi) == group1A)
  row.choice.1B <- which(rownames(autosome.genome.rolling.stats$pi) == group1B)
  row.choice.2A <- which(rownames(autosome.genome.rolling.stats$pi) == group2A)
  row.choice.2B <- which(rownames(autosome.genome.rolling.stats$pi) == group2B)
  vector1A <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1A,])
  vector1B <- as.vector(autosome.genome.rolling.stats$pi[row.choice.1B,])
  mean_pi_1 <- (vector1A + vector1B) / 2
  vector2A <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2A,])
  vector2B <- as.vector(autosome.genome.rolling.stats$pi[row.choice.2B,])
  mean_pi_2 <- (vector2A + vector2B) / 2
  Dxy.row.choice_1 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_1)
  Dxy.row.choice_2 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_2)
  Dxy_1 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_1,])
  Dxy_2 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_2,])
  mean_pi_over_Dxy_1 <- mean_pi_1 / Dxy_1
  mean_pi_over_Dxy_2 <- mean_pi_2 / Dxy_2
  quartz(title=paste0("Scatterplot of mean windowed pi/Dxy, based on compiled chromosomes, between", group1A,"_", group1B, " vs ", group2A,"_", group2B, sep=""), width=4.25, height=5)
  plot(x=mean_pi_over_Dxy_1, y=mean_pi_over_Dxy_2, 
       xaxp = c(0, 0.5, 5), yaxp=c(0, 0.5, 5),
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mar=c(3,3,0,0), mgp=c(3,0.5,0))
  # add labels:
  label.size <- 1.2
  title(xlab=bquote("Mean " * italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(group1A)) * " to " * italic(.(group1B)))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom("Mean " * italic(pi) / italic(D)[xy] * ", ") * italic(.(group1A)) * " to " * italic(.(group1B))), line=2.5, cex.lab=label.size, col.lab=color1)
  title(ylab=bquote("Mean " * italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(group2A)) * " to " * italic(.(group2B)))), line=2, cex.lab=label.size)
  title(ylab=bquote(phantom("Mean " * italic(pi) / italic(D)[xy] * ", ") * italic(.(group2A)) * " to " * italic(.(group2B))), line=2, cex.lab=label.size, col.lab=color2)
  test <- cor.test(mean_pi_over_Dxy_1, mean_pi_over_Dxy_2)
  return(test)
}


# Scatterplot of Dxy vs. mean_pi using compiled genome-wide info
# Color the points based on Fst
# Very cool graph
plotDxy_MeanPi <- function(group1, group2, cor.method,
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
  color_scale <- round(WC84_Fst.vector * 100) # / max.WC84_Fst)
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
  title(xlab=expression("Between-group differentiation (" * italic('D')[xy] * ")"), line=2.5, cex.lab=label.size)
  title(ylab=expression("Mean within-group variation (" * italic(pi) * ")"), line=2, cex.lab=label.size)
  test <- cor.test(Dxy, mean_pi, method=cor.method)
  return(test)
}


# Make a scatterplot of Dxy/pi vs. WC84_Fst:
plotDxyOverMeanPi_Fst <- function(group1, group2,
                           autosome.genome.rolling.stats) {
  groups.for.graph <- paste0(group1, "_", group2)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  mean_pi <- (autosome.genome.rolling.stats$pi[row.choice.1,] + autosome.genome.rolling.stats$pi[row.choice.2,]) / 2
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
  Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  WC84_Fst.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == groups.for.graph)
  WC84_Fst.vector <- autosome.genome.rolling.stats$WC84_Fst[WC84_Fst.row.choice,]
  # color points as a gradient according to WC84_Fst:
  color.palette <- colorRampPalette(c("blue", "red"))(100)  # makes a color scale with 100 colors
  max.WC84_Fst <- max(is.finite(WC84_Fst.vector))  # the next line will cause the max.WC84_Fst to be set to the top color
  color_scale <- round(WC84_Fst.vector * 100 / max.WC84_Fst)
  color_scale[color_scale<1] <- 1
  colors <- color.palette[color_scale]  # makes a list of colors according to Fst values of the windows
  quartz(title=paste0("Scatterplot of windowed Dxy/(mean_pi) vs. WC84_Fst between ", group1,"_", group2, sep=""), width=6, height=7.5)
  plot(Dxy/mean_pi, WC84_Fst.vector, col=colors)
}


# Make a scatterplot of Dxy vs. pi of one group:
plotDxy_Pi <- function(group1, group2, pi_group, cor.method,
                       autosome.genome.rolling.stats) {
  groups.for.graph <- paste0(group1, "_", group2)
  pi.row.choice <- which(rownames(autosome.genome.rolling.stats$pi) == pi_group)
  pi <- autosome.genome.rolling.stats$pi[pi.row.choice,]
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
  Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  WC84_Fst.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == groups.for.graph)
  WC84_Fst.vector <- autosome.genome.rolling.stats$WC84_Fst[WC84_Fst.row.choice,]
  # color points as a gradient according to WC84_Fst:
  color.palette <- colorRampPalette(c("blue", "red"))(100)  # makes a color scale with 100 colors
  max.WC84_Fst <- max(is.finite(WC84_Fst.vector))  # the next line will cause the max.WC84_Fst to be set to the top color
  color_scale <- round(WC84_Fst.vector * 100 / max.WC84_Fst)
  color_scale[color_scale<1] <- 1
  colors <- color.palette[color_scale]  # makes a list of colors according to Fst values of the windows
  quartz(title=paste0("Scatterplot of windowed Dxy for ", group1,"_", group2, " vs. pi in ", pi_group, sep=""), width=6, height=6)
  plot(Dxy, pi, col=colors, xlim=c(0, max(Dxy)*1.04), ylim=c(0, max(pi)*1.04), xlab=NA, ylab=NA)
  lines(c(0,1), c(0,1))
  label.size <- 1.1
  title(xlab=bquote("Between-group differentiation (" * italic('D')[xy] * ")"), line=3, cex.lab=label.size)
  title(ylab=bquote(italic(pi) * " in " * italic(.(pi_group))), line=2.5, cex.lab=label.size)
  test <- cor.test(Dxy, pi, method=cor.method)
  return(test)
}


# Scatterplot of Pi/Dxy in/between two groups
# points colored by Fst (red indicates higher Fst)
plotPiOverDxy_corr <- function(group1, group2, cor.method,
                               autosome.genome.rolling.stats) {
  groups.for.graph <- paste0(group1, "_", group2)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  mean_pi <- (autosome.genome.rolling.stats$pi[row.choice.1,] + autosome.genome.rolling.stats$pi[row.choice.2,]) / 2
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
  Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  WC84_Fst.row.choice <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == groups.for.graph)
  WC84_Fst.vector <- autosome.genome.rolling.stats$WC84_Fst[WC84_Fst.row.choice,]
  # color points as a gradient according to WC84_Fst:
  color.palette <- colorRampPalette(c("blue", "red"))(100)  # makes a color scale with 100 colors
  max.WC84_Fst <- max(is.finite(WC84_Fst.vector))  # the next line will cause the max.WC84_Fst to be set to the top color
  color_scale <- round(WC84_Fst.vector * 100 / max.WC84_Fst)
  color_scale[color_scale<1] <- 1
  colors <- color.palette[color_scale]  # makes a list of colors according to Fst values of the windows
  pi_1_over_Dxy <- pi_1 / Dxy
  pi_2_over_Dxy <- pi_2 / Dxy
  quartz(title=paste0("Scatterplot of pi_1/Dxy vs. pi_2/Dxy for both ", group1,"_", group2, sep=""), width=6, height=6)
  plot(pi_1_over_Dxy, pi_2_over_Dxy, col=colors, xlim=c(0, max(pi_1_over_Dxy)*1.04), ylim=c(0, max(pi_2_over_Dxy)*1.04))
  test <- cor.test(pi_1_over_Dxy, pi_2_over_Dxy, method=cor.method)
  return(test)
}


# Graph of 3-way comparison of windowed pi / Dxy_max
PiOverDxyMaxMatrixPlot3 <- function(group1, group2, group3, cor.method, 
                                    autosome.genome.rolling.stats)  {
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$pi) == group3)
  pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  pi_3 <- autosome.genome.rolling.stats$pi[row.choice.3,]
  group.pair_12 <- paste0(group1, "_", group2)
  group.pair_13 <- paste0(group1, "_", group3)
  group.pair_23 <- paste0(group2, "_", group3)
  Dxy.row.choice_12 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_12)
  Dxy.row.choice_13 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_13)
  Dxy.row.choice_23 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_23)
  Dxy_12 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_12,])
  Dxy_13 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_13,])
  Dxy_23 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_23,])
  Dxy_max <- pmax(Dxy_12, Dxy_13, Dxy_23)
  pi_1_over_Dxy_max <- pi_1 / Dxy_max
  pi_2_over_Dxy_max <- pi_2 / Dxy_max
  pi_3_over_Dxy_max <- pi_3 / Dxy_max
  quartz(title=paste0("Scatterplot matrix of windowed pi / Dxy_max, based on compiled chromosomes, ", group1," vs ", group2, " vs ", group3, sep=""), width=6, height=6)
  pairs(cbind(pi_1_over_Dxy_max, pi_2_over_Dxy_max, pi_3_over_Dxy_max), labels=c(group1, group2, group3), pch=16, cex=0.25,
        lower.panel=NULL, diag.panel=NULL)
  test1 <- cor.test(pi_1_over_Dxy_max, pi_2_over_Dxy_max, method=cor.method)
  test2 <- cor.test(pi_1_over_Dxy_max, pi_3_over_Dxy_max, method=cor.method)
  test3 <- cor.test(pi_2_over_Dxy_max, pi_3_over_Dxy_max, method=cor.method) 
  return(list(test1=test1, test2=test2, test3=test3))
}


# Make 3 scatterplots of pi/Dxy_max in a line:
plots3PiOverDxyMax_corr <- function(groups, colors, cor.method,
                                    autosome.genome.rolling.stats) {
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == groups[1])
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == groups[2])
  row.choice.3 <- which(rownames(autosome.genome.rolling.stats$pi) == groups[3])
  pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  pi_3 <- autosome.genome.rolling.stats$pi[row.choice.3,]
  group.pair_12 <- paste0(groups[1], "_", groups[2])
  group.pair_13 <- paste0(groups[1], "_", groups[3])
  group.pair_23 <- paste0(groups[2], "_", groups[3])
  Dxy.row.choice_12 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_12)
  Dxy.row.choice_13 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_13)
  Dxy.row.choice_23 <- which(rownames(autosome.genome.rolling.stats$Dxy) == group.pair_23)
  Dxy_12 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_12,])
  Dxy_13 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_13,])
  Dxy_23 <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice_23,])
  Dxy_max <- pmax(Dxy_12, Dxy_13, Dxy_23)
  vector1 <- pi_1 / Dxy_max
  vector2 <- pi_2 / Dxy_max
  vector3 <- pi_3 / Dxy_max
  quartz(title=paste0("Scatterplots of pi/Dxy_max in three comparisons"), width=8, height=3)
  par(oma=c(1,3,1,3))  # set outer margins
  par(mfrow=c(1,3))
  label.size <- 1.5
  high.lim <- round(max(c(vector1, vector2, vector3)), 2) + 0.01
  # first plot
  plot(x=vector1, y=vector2, xlim=c(0,high.lim), xaxp=c(0,1,2), ylim=c(0,high.lim), yaxp=c(0,1,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(groups[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic(pi) / italic(D)[xy] * ", ") * italic(.(groups[1]))), line=2.5, cex.lab=label.size, col.lab=colors[1])
  title(ylab=bquote(italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(groups[2])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi) / italic(D)[xy] * ", ") * italic(.(groups[2]))), line=2.5, cex.lab=label.size, col.lab=colors[2])
  # second plot
  plot(x=vector1, y=vector3, xlim=c(0,high.lim), xaxp=c(0,1,2), ylim=c(0,high.lim), yaxp=c(0,1,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(groups[1])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic(pi) / italic(D)[xy] * ", ") * italic(.(groups[1]))), line=2.5, cex.lab=label.size, col.lab=colors[1])
  title(ylab=bquote(italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(groups[3])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi) / italic(D)[xy] * ", ") * italic(.(groups[3]))), line=2.5, cex.lab=label.size, col.lab=colors[3])
  # third plot
  plot(x=vector2, y=vector3, xlim=c(0,high.lim), xaxp=c(0,1,2), ylim=c(0,high.lim), yaxp=c(0,1,2), 
       asp=1, pch=16, cex=0.25, cex.axis=0.9, tcl=-0.5, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
  title(xlab=bquote(italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(groups[2])))), line=2.5, cex.lab=label.size)
  title(xlab=bquote(phantom(italic(pi) / italic(D)[xy] * ", ") * italic(.(groups[2]))), line=2.5, cex.lab=label.size, col.lab=colors[2])
  title(ylab=bquote(italic(pi) / italic(D)[xy] * ", " * phantom(italic(.(groups[3])))), line=2.5, cex.lab=label.size)
  title(ylab=bquote(phantom(italic(pi) / italic(D)[xy] * ", ") * italic(.(groups[3]))), line=2.5, cex.lab=label.size, col.lab=colors[3])
  test1 <- cor.test(vector1, vector2, method=cor.method)
  test2 <- cor.test(vector1, vector3, method=cor.method)
  test3 <- cor.test(vector2, vector3, method=cor.method) 
  return(list(test1=test1, test2=test2, test3=test3))
}


# Compile variant numbers ----
# Reads files containing per-bp info for each chromosome;
# requires that "Fst_among" is in the WC84_Fst matrix;
# returns a summary of variant and invariant bp per chromosome
compileVariantNumbers <- function(base.name, tag.name, chromosomes.to.summarize) {
  loci.in.chr <- NULL
  snps.in.chr <- NULL
  invariant.loci.in.chr <- NULL
  for (i in 1:length(chromosomes.to.summarize)) {
    chr <- chromosomes.to.summarize[i] 
    region.text <- paste0("Chr",chr,"_whole")
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print(paste0("Loaded saved summary site stats for ",region.text, sep=""))
    loci.in.chr[i] <- length(WC84_Fst[1,])
    snps.in.chr[i] <- length(which(!is.na(WC84_Fst[which(rownames(WC84_Fst)=="Fst_among"),])))
    invariant.loci.in.chr[i] <- length(which(is.na(WC84_Fst[which(rownames(WC84_Fst)=="Fst_among"),])))
  }
  num_loci_summary <- data.frame(loci=loci.in.chr, SNPs=snps.in.chr, invariant=invariant.loci.in.chr, row.names=chromosomes.to.summarize)
  return(num_loci_summary)
}


# MeanPi_Dxy autosome vs. Z plot ----
# uses compiled data produced by function "compileWindowedStats";
# for both autsome and Z
plotMeanPi_Dxy.autosome_Z <- function(group1, group2,
                                      autosome.genome.rolling.stats, Z.rolling.stats) {
  groups.for.graph <- paste0(group1, "_", group2)
  row.choice.1 <- which(rownames(autosome.genome.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(autosome.genome.rolling.stats$pi) == group2)
  autosome.pi_1 <- autosome.genome.rolling.stats$pi[row.choice.1,]
  autosome.pi_2 <- autosome.genome.rolling.stats$pi[row.choice.2,]
  autosome.mean_pi <- (autosome.pi_1 + autosome.pi_2) / 2
  Dxy.row.choice <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
  autosome.Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice,])
  # Z calcs:
  row.choice.1 <- which(rownames(Z.rolling.stats$pi) == group1)
  row.choice.2 <- which(rownames(Z.rolling.stats$pi) == group2)
  Z.pi_1 <- Z.rolling.stats$pi[row.choice.1,]
  Z.pi_2 <- Z.rolling.stats$pi[row.choice.2,]
  Z.mean_pi <- (Z.pi_1 + Z.pi_2) / 2
  Dxy.row.choice <- which(rownames(Z.rolling.stats$Dxy) == groups.for.graph)
  Z.Dxy <- as.vector(Z.rolling.stats$Dxy[Dxy.row.choice,])
  # make the figure:
  quartz(title=paste0("Autosomes vs. Z: Scatterplot of windowed mean pi vs. Dxy between ", group1,"_", group2, sep=""), width=6, height=6)
  par(oma=c(3,3,1,1))  # set outer margins
  zones <- matrix(c(4,0,0,2,0,0,1,3,5), ncol=3, byrow=TRUE)  # numbers in matrix give order of plotting
  layout(zones, widths=c(3/5,1/5,1/5), heights=c(1/5,1/5,3/5))
  xlimits <- c(0, round(max(c(autosome.Dxy,Z.Dxy))*1.1, digits=3))  # rounds to 3 decimal places
  ylimits <- c(0, round(max(c(autosome.mean_pi,Z.mean_pi))*1.1, digits=3))
  limits <- c(min(c(xlimits,ylimits)), max(xlimits,ylimits)) # this to make axes have same limits
  xhist <- hist(autosome.Dxy, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
  chrZ.xhist <- hist(Z.Dxy, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
  yhist <- hist(autosome.mean_pi, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
  chrZ.yhist <- hist(Z.mean_pi, breaks=seq(limits[1], limits[2], by=0.0005), plot=FALSE)
  top <- max(c(xhist$counts, yhist$counts))
  chrZ.top <- max(c(chrZ.xhist$counts, chrZ.yhist$counts))
  par(mar=c(3,3,1,1))  # specifies number of lines around plot (bottom, left, top right)
  plot(autosome.Dxy,autosome.mean_pi, col="grey70", xlim=limits, 
       ylim=limits, pch=16, cex=0.3, asp=1)
  lines(c(0,1), c(0,1))
  points(Z.Dxy, Z.mean_pi, col='blue', pch=16, cex=0.5, asp=1)
  par(mar=c(0,3,0,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), col="grey70", space=0)
  axis(side=2, at=c(0,200))
  par(mar=c(3,0,1,0.5))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), col="grey70", space=0, horiz=TRUE)
  axis(side=1, at=c(0,200))
  par(mar=c(0,3,1,1))
  barplot(chrZ.xhist$counts, axes=FALSE, ylim=c(0, top), space=0, col="blue")
  axis(side=2, at=c(0,200))
  par(mar=c(3,0,1,1))
  barplot(chrZ.yhist$counts, axes=FALSE, xlim=c(0, top), space=0, col="blue", horiz=TRUE)
  axis(side=1, at=c(0,200))
  par(oma=c(3,3,0,0))
  label.size <- 1.2
  mtext(expression("Between-group differentiation (" * italic('D')[xy] * ")"), side=1, line=0.5, outer=TRUE, 
        at=1.5/5)
  mtext(expression("Mean within-group variation (" * italic(pi) * ")"), side=2, line=0, outer=TRUE, 
        at=1.5/5)
  mtext("Windows", side=1, line=-0.5, outer=TRUE, cex=0.8, 
        at=3.15/5)
  mtext("Windows", side=1, line=-0.5, outer=TRUE, cex=0.8, 
        at=4.1/5)
  mtext("Windows", side=2, line=-0.5, outer=TRUE, cex=0.8, 
        at=3.15/5)
  mtext("Windows", side=2, line=-0.5, outer=TRUE, cex=0.8, 
        at=4.1/5)
  # t-tests:
  test.mean_pi <- t.test(autosome.mean_pi, Z.mean_pi)
  test.Dxy <- t.test(autosome.Dxy, Z.Dxy)
  return(list(test.mean_pi=test.mean_pi, test.Dxy=test.Dxy))
}
