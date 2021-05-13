#This script was used to find the mean and median values of pi_b, and pi_w per chromosome.
#And calculate averages of the windowed Fst per chromosome
#This data was used to produce the plots of Fst, pi_b, and pi_w vs chromosome length

#setup
setwd("~/Desktop/Wrens/pawr_wiwr_genomics")
source ("~/Desktop/Wrens/pawr_wiwr_genomics/PAWR_WIWR_genomics_Rproject/Genomics_R_Functions_changeddxytopib.R")


##### Calculate per-chromosome pi_w and pi_b #####

#select chromosomes to analyze (all of them)
chromosomes.to.analyze <- c("1", "1A", "2", "3","4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15","17","18","19","20","21","22","23","24","25","26","27","28","Z")

#more setup for the loop. We are loading pre-calculated data from previous runs to save time
Analysis_set <- 2  # 1: all samples, only SNPs;    2: all samples, with invariant sites

if (Analysis_set == 1) {
  # nothing
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
  group.colors.WC84_Fst <- c("purple", "blue", "red", "grey", "orange")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("PAWR", "WIWR", "Hybrid", "MAWR")
  group.colors.pi <- c("blue", "red", "purple", "grey")
  
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

#initialize empty vectors to hold data
meanFst <- NA
medianFst <- NA
meanpib <- NA
meanpiPAWR <- NA
meanpiWIWR <- NA

# MAIN LOOP ----
# -----
#I have deleted much of the main loop to save time, I only want to get per-chromosome means and medians
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
  
  #find out which row contains the data of interest
  groups.to.compare <- "PAWR_WIWR"
  row.choice <- which(rownames(WC84_Fst) == groups.to.compare)
  #could select an fst cutoff if desired
  #WC84_Fst.cutoff <- 0.9
  #selection <- (WC84_Fst[row.choice,] > WC84_Fst.cutoff) & (!is.na(WC84_Fst[row.choice,])) & (pos$position > start.pos) & (pos$position < end.pos)
  print(chr)
  #print(sum(selection))
  #print(mean(WC84_Fst[row.choice,], na.rm=TRUE)) #this would give you the mean of individual sites, but we are instead using the W&C method to get windowed estimates of Fst
  #print(median(WC84_Fst[row.choice,], na.rm=TRUE))
  meanFst[i]<-print(mean(rolling.mean.WC84_Fst[row.choice,], na.rm=TRUE))
  medianFst[i]<-print(median(rolling.mean.WC84_Fst[row.choice,], na.rm=TRUE))
  meanpib[i]<-print(mean(Dxy[row.choice,], na.rm=TRUE))
  #print(median(Dxy[row.choice,], na.rm=TRUE)) #median Dxy is zero, not informative
  meanpiPAWR[i]<-print(mean(site_pi_nb["PAWR",], na.rm=TRUE))
  #print(median(site_pi_nb["PAWR",], na.rm=TRUE)) #median pi is zero, not informative
  meanpiWIWR[i]<-print(mean(site_pi_nb["WIWR",], na.rm=TRUE))
  #print(median(site_pi_nb["WIWR",], na.rm=TRUE)) #median pi is zero, not informative
  
}   # End main loop ----

#view results. These could be loaded as vectors in the future to save time recalculating them.
dput(chromosomes.to.analyze)
# c("1", "1A", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", 
#   "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", 
#   "23", "24", "25", "26", "27", "28", "Z")

dput(meanFst)
#c(0.229475273314787, 0.250764670471225, 0.244278441233746, 0.25725709367356, 
# 0.235117506101421, 0.267735386972732, 0.233919773854388, 0.245922848363538, 
# 0.260965021841355, 0.266642488455805, 0.260506336749497, 0.247005219197547, 
# 0.251037632398637, 0.300423324110158, 0.215662563160699, 0.26932033823055, 
# 0.260181757827567, 0.259218259429947, 0.287076012988077, 0.240194611606385, 
# 0.239625746455673, 0.214767565233846, 0.295557094751859, 0.188897144889735, 
# 0.273534641946205, 0.420220541497801, 0.227513359497241, 0.383995361436798, 
# 0.335600628865456, 0.359387082320617)

dput(medianFst)
# c(0.214956176588721, 0.23302674409908, 0.215521643300931, 0.222234084048413, 
#   0.218881418702786, 0.203203111495488, 0.205000280060415, 0.20327891364625, 
#   0.23902219522425, 0.237904979648208, 0.210922122115021, 0.193513547910174, 
#   0.21486144987219, 0.229144306644943, 0.213038813483801, 0.234212274767008, 
#   0.232814335213418, 0.206731516567913, 0.229603593442156, 0.219208754221705, 
#   0.199642973178636, 0.188595386027539, 0.269750621074986, 0.154936622560533, 
#   0.230036943850693, 0.403646396214235, 0.246754708040892, 0.398905347251677, 
#   0.299617073417575, 0.349786347415645)

dput(meanpib)
# c(0.0040780920000196, 0.00400438764269604, 0.00413041287880188, 
#   0.00407750504033084, 0.00429125072217047, 0.00358987479247667, 
#   0.00420358284021871, 0.0038392728963954, 0.00417831466441949, 
#   0.00423461042141485, 0.00395454198946797, 0.00364450084295963, 
#   0.00385682088227635, 0.00410697605341887, 0.00347078435468122, 
#   0.00317086244357316, 0.00359210731457206, 0.00340105421118581, 
#   0.00371815558284978, 0.00375539527419549, 0.00369645000323785, 
#   0.00346531588695395, 0.00372621143853866, 0.0032141778528496, 
#   0.00344278580033459, 0.00280673853872536, 0.0036999729909087, 
#   0.00328526097341095, 0.00337023889823053, 0.00395475131055895
# )

dput(meanpiPAWR)
# c(0.00285245833260542, 0.00266787929181686, 0.0028095793864253, 
#   0.00274180896159134, 0.00299907266661534, 0.00239449542565647, 
#   0.00289832496194629, 0.00264206197566118, 0.0026867206916576, 
#   0.00280340429698883, 0.00273692752162319, 0.00243466361243029, 
#   0.00262365839408402, 0.00252128628574022, 0.00241491504280418, 
#   0.00201154758179498, 0.00240907264903217, 0.00230141034785786, 
#   0.00245084443165617, 0.00245734991919961, 0.00249660323164725, 
#   0.00245141641460447, 0.00209017092214821, 0.00237861518339936, 
#   0.00213346097471062, 0.00115095286888068, 0.00272408413555423, 
#   0.00182080769297508, 0.00183515336670245, 0.00212554773470473
# )

dput(meanpiWIWR)
# c(0.00356529469913947, 0.00339694021783188, 0.00358385166741606, 
#   0.00338046069363756, 0.0036668161879359, 0.00308377394231129, 
#   0.00367452873711146, 0.00325902803297958, 0.00355240839492167, 
#   0.00355986642656444, 0.00320852205019791, 0.0031252796035144, 
#   0.0032165498266832, 0.00334498332979905, 0.00305598728048008, 
#   0.00256966076143401, 0.00302341516553586, 0.00289007842630943, 
#   0.00287441684792372, 0.00344996315038529, 0.0032197449291388, 
#   0.00308077303209173, 0.00352653651211573, 0.00288571843510117, 
#   0.00280671046986167, 0.00197836228082518, 0.00299180495471171, 
#   0.0022658751043063, 0.00237665723118792, 0.00303965028271846)


##### Calculate per-chromosome Fst #####
#Next I need to calculate Fst per chromosome. We could use the meanFst vector that we just calculated, but this is not ideal for two reasons
#1) The Weir&Cockerham 1984 method for calculating Fst with multiple sites does not take the mean of the sites, it takes the mean of the numerator divided by the mean of the denominator, which gives a slightly different value
#2) We calculated Fst in windows containing 10,000 genotyped positions, but SNP density varies among regions, so that each rolling mean window actually is calculated with different numbers of SNPs. If we took the mean of the windows with the invariants included, each SNP would not be weighted equally.

#do not need invariants for Fst calculation
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
  # ID	location	group	Fst_group	plot_order
  metadata.file <- "PAWRWIWR_Metadata.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 77  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  #window_size <- 1000
  window_size <- 10000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR")
  group.colors <- c("blue", "red")
  group_count <- length(groups)
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR")
  group.colors.WC84_Fst <- c("purple")
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
  # ID  location  group Fst_group plot_order
  metadata.file <- "PAWRWIWR_Metadata.txt"  # Note the first line in this file corrects an error: AD29A07 should be ED29A07
  num.individuals <- 77  # specify number of individuals in file (good for error checking)
  # specify window size (number of bp with info) and step size
  window_size <- 10000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("PAWR", "WIWR")
  group.colors <- c("blue", "red")
  group_count <- length(groups)
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("PAWR_WIWR")
  group.colors.WC84_Fst <- c("purple")
}


# Option to focus on a region of chromosome ----
focus.region <- F  # choose T for a subset of the chromosome, F for the whole thing)

calculate_or_load_stats <- 2  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) load per-site and windowed data from file

# Option to save the per-site stats
saveSiteInfo <- F # TRUE     # TRUE   # If TRUE, will save a file for per-site stats
saveWindowedStats <- F # TRUE     #TRUE   # If TRUE, will save a file for per-window stats
load.rolling.means <- F   #FALSE     #   FALSE    # If TRUE, will load rolling mean data (rather than calculate)
locations <- read.table(paste0(metadata.file), sep = "\t", header=TRUE)
num_loc_cols <- length(locations[1,])
meanChrFst <- NA

# calculate chromosome-wide Fst according to according to Weir&Cockerham1984 (with sample size and pop number correction),
# calculated as windowed numerator over windowed denominator,
getChrWC84_Fst <- function(WC84_Fst_numerator, WC84_Fst_denominator, row.choice){
  mean.WC84_Fst_numerator <- mean(WC84_Fst_numerator[row.choice,], na.rm = TRUE)
  mean.WC84_Fst_denominator <- mean(WC84_Fst_denominator[row.choice,], na.rm = TRUE)
  mean.WC84_Fst <- mean.WC84_Fst_numerator / mean.WC84_Fst_denominator
  return(mean.WC84_Fst)
}

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
    
  } else if (calculate_or_load_stats==2 | calculate_or_load_stats==3) {
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print("Loaded saved summary stats")
  }
  
  groups.to.compare <- "PAWR_WIWR"
  row.choice <- which(rownames(WC84_Fst) == groups.to.compare)
  chr_Fst<- getChrWC84_Fst(WC84_Fst_numerator, WC84_Fst_denominator, row.choice)
  print(paste("Fst for ", chr, ":", chr_Fst))
  meanChrFst[i] <- chr_Fst
}   
# End main loop ----
#get Fst results
dput(meanChrFst)
#c(0.26732052846208, 0.282452937077198, 0.276737944986152, 0.298326392592223, 
#0.271534255357431, 0.322040678657465, 0.274527465505898, 0.283663901294031, 
#0.307231441167758, 0.308464995011577, 0.32448503904949, 0.295691742405375, 
#0.29519094771614, 0.347426634293701, 0.264496987979719, 0.322741772723118, 
#0.323971159361559, 0.296782557484343, 0.336208211586033, 0.295453420308393, 
#0.28563932980528, 0.275672105817445, 0.470743945304538, 0.266607219396253, 
#0.346334635553296, 0.671081035421665, 0.253772004669021, 0.466203091530677, 
#0.446884728034056, 0.417505579741889)

#just to check that the window size datasets are identical, I tried again with the " window_size <- 10000" dataset and the values are completely identical as expected
#c(0.26732052846208, 0.282452937077198, 0.276737944986152, 0.298326392592223, 
#0.271534255357431, 0.322040678657465, 0.274527465505898, 0.283663901294031, 
#0.307231441167758, 0.308464995011577, 0.32448503904949, 0.295691742405375, 
#0.29519094771614, 0.347426634293701, 0.264496987979719, 0.322741772723118, 
#0.323971159361559, 0.296782557484343, 0.336208211586033, 0.295453420308393, 
#0.28563932980528, 0.275672105817445, 0.470743945304538, 0.266607219396253, 
#0.346334635553296, 0.671081035421665, 0.253772004669021, 0.466203091530677, 
#0.446884728034056, 0.417505579741889)

##### Plot results #####
#make a vector containing the chromosome lengths, in the same order that they are listed in the chromosomes.to.analyze vector
ChrSizes <- c(119.8, 74.8, 157.4, 115.7, 70.3, 21.2, 64.6, 37.2, 39.3, 32, 26.8, 21.3, 21.7, 21.9, 18.6, 17.4, 14.9, 12.4, 13.1, 11.9, 15.6, 8.1, 5.7, 7.9, 8, 1.3, 4.9, 4.6, 5, 74.6)

#plot chromosome length vs pi
quartz()
plot(ChrSizes, meanpib, col="black", main = " ",
     xlab = "Chromosome Length (Mb)", ylab = "Nucleotide Diversity",
     pch = 19, frame = FALSE, xlim=c(0, 120), ylim=c(0, 0.005))
#add pi to the plot
points(ChrSizes, meanpiPAWR, col="blue", pch=19)
points(ChrSizes, meanpiWIWR, col="red", pch=19)
#highlight the points that belong to the Z chromosome by plotting over them in gold (the values for the Z chromosome can be found as the last entry of each vector)
points(74.6, 0.002125548, col="orange", pch=19)
points(74.6, 0.00303965, col="orange", pch=19)
points(74.6, 0.003954751, col="orange", pch=19)
#add rings of colour around the gold points so the reader can tell which statistic they belong to
points(74.6, 0.002125548, col="blue")
points(74.6, 0.00303965, col="red")
points(74.6, 0.003954751, col="black")

#plot chromosome length vs fst using the means of 10kb rolling means (I want to compare)
quartz()
plot(ChrSizes, meanFst, col="black", main = " ",
     xlab = "Chromosome Length (Mb)", ylab = "FST",
     pch = 19, frame = FALSE, xlim=c(0, 120), ylim=c(0, 0.8))
#highlight the points that belong to the Z chromosome by plotting over them in gold (the values for the Z chromosome can be found as the last entry of each vector)
points(74.6, 0.359387082320617, col="orange", pch=19)
#add rings of colour around the gold points so the reader can tell which statistic they belong to
points(74.6, 0.359387082320617, col="black")

#plot chromosome length vs fst using the windowed estimate for the entire chromosome as a single window (this is what I will publish)
quartz()
plot(ChrSizes, meanChrFst, col="black", main = " ",
     xlab = "Chromosome Length (Mb)", ylab = "FST",
     pch = 19, frame = FALSE, xlim=c(0, 120), ylim=c(0, 0.8))
#highlight the points that belong to the Z chromosome by plotting over them in gold (the values for the Z chromosome can be found as the last entry of each vector)
points(74.6, 0.417505579741889, col="orange", pch=19)
#add rings of colour around the gold points so the reader can tell which statistic they belong to
points(74.6, 0.417505579741889, col="black")

#I have saved the two raw figure panels into the Github folder 5_GenomeScans_Plots/means_per_chromosome
#I will combine the panels into a single figure in Illustrator, and try different fonts for the publication
