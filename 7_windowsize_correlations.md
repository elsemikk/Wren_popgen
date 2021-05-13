It was suggested during review to check whether the statistics (FST, pi_w, or pi_b) are correlated with window size in this dataset. Here is the code I used to check these correlations. This code requires the data to be read into R in the previous step (outlined in the GenomeScans page).
```R
#First, find the sizes of each of the windows
getWindowSize <- function(pos, window_size, step_size){
  rolling.min.pos <- rollapply(pos$position, width=window_size, FUN=min, by=step_size, align="left")
  rolling.max.pos <- rollapply(pos$position, width=window_size, FUN=max, by=step_size, align="left")
  rolling.size <- rolling.max.pos - rolling.min.pos
  return(rolling.size)
}
rolling.size <- getWindowSize(pos, window_size, step_size)

#after loading the data (from the genome scan analysis), get vectors to hold Fst, dxy, and pi
mean_WC84_Fst=rolling.mean.WC84_Fst[1,]
mean_Dxy=rolling.mean.Dxy[1,]
mean_pi=rolling.mean.pi[1,]
#plot the relationship between Fst, pi_b (Dxy), and pi_w with window size
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

#Here is how to loop through each chromosome to make these plots for each one:
window_size <- 10000
step_size <- window_size  # could change if wanting overlapping windows
chromosomes.to.analyze <- c("1", "1A", "2", "3","4","4A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15","17","18","19","20","21","22","23","24","25","26","27","28","Z")

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
  quartz()
  print(ggplot(as.data.frame(mean_WC84_Fst), aes(x=rolling.size, y=mean_WC84_Fst))+
    geom_point()+
    theme_classic()+
    xlab("Window Size")+
    ylab("Fst")+
    stat_smooth(method=lm, aes(x=rolling.size, y=mean_WC84_Fst)))

  quartz()
  print(ggplot(as.data.frame(mean_pi), aes(x=rolling.size, y=mean_pi))+
    geom_point()+
    theme_classic()+
    xlab("Window Size")+
    ylab("pi_w")+
    stat_smooth(method=lm, aes(x=rolling.size, y=mean_pi)))

  quartz()
  print(ggplot(as.data.frame(mean_Dxy), aes(x=rolling.size, y=mean_Dxy))+
    geom_point()+
    theme_classic()+
    xlab("Window Size")+
    ylab("pi_b")+
    stat_smooth(method=lm, aes(x=rolling.size, y=mean_Dxy)))
} 
```

I have placed the resulting plots into this repository, into the folder `5_GenomeScans_Plots/Windowsize_correlations`.