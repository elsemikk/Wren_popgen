This script provides the calculation of heterozygosity at the called genotypes of each sample, as well as the heterozygosity of each sample at the Z chromosome, which was used to identify the sex of each sample (as males have two copies of the Z, while females have only one Z and should show very low heterozygosity at Z-linked loci except for genotyping/mapping errors and sites within the pseudoautosomal region where the Z and W chromosomes undergo recombination)

This script was used to estimate autosomal heterozygosity (ie. excluding chromosome Z):
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
```


In order to identify males and females, we can compare heterozygosities on the Z chromosome. 

First, I loaded the data as per the `5_GenomeScans` pipeline, loading the whole dataset (rather than `chromosomes.to.analyze <- c("autoinfSNP")`). Since chromsome Z was listed last in the list of chromosomes, the data for chromosome Z was left in the `geno` variable (with 422,161 sites). I could then run the heterozygosity calculation as above, and make a plot. This plot is provided in the `5_GenomeScans_Plots` folder as `chrZ_heterozygosity.pdf`.
```R
#set up an empty vector to hold the data
heterozygosities <- NULL
#For each row of genotypes, calculate heterozygosity (the number of "1" genotype calls divided by the total number of non-missing calls for that sample))
for(i in 1:nrow(geno)){
  heterozygosities$heterozygosity[i]<- sum(geno[i,]=="1", na.rm=T)/sum(!is.na(geno[i,]))
}
heterozygosities$sample <- ind
heterozygosities$location <- locations
heterozygosities <- as.data.frame(heterozygosities)

#make another column listing the sample's sex. Conveniently, all females are below 0.008 and the males are above 0.008
heterozygosities$sex <- "Male"
heterozygosities$sex[heterozygosities$heterozygosity<0.0008] <- "Female"

#plot the data
library(dplyr)
heterozygosities %>%
filter(location.group != "low_read_WIWR") %>% #remove the failed sample that has (virtually) no data
  ggplot(., aes(x=sex, y=heterozygosity)) + #create a ggplot
  theme_classic() + #remove background colour and extra lines
  theme(axis.text.y = element_text(family="Times", face="plain", colour="black", size=10), axis.text.x = element_text(family="Times", face="plain", colour="black", size=10), axis.title.x = element_text(family="Times", face="plain", colour="black", size=14), axis.title.y = element_text(family="Times", face="plain", colour="black", size=14), legend.text = element_text(family="Times", colour="black", size=10), legend.title = element_text(family="Times", face="plain", colour="black", size=12))+ #define the text fonts and sizes for presentation purposes
  xlab("Sex")+ #write x axis label
  ylab("Z Chromosome Heterozygosity")+ #write y axis label
  labs(colour="Species")+ #make legend label
  scale_color_manual(values=c("purple", "green", "blue", "red"), labels = c("Hybrid", "Marsh", "Pacific", "Winter")) + #set colours to match population colour scheme
  geom_boxplot(outlier.size=0) + #I am removing the outlier dots because they are redundant with the jitter plot and can be confused for additional data points
  geom_jitter(binaxis='y', stackdir='center', dotsize=1, aes(color = location.group)) #plot the data points on top of the boxplot

```
Here is the result of the calculation of heterozygosities. Note that aong males, the  WIWR have higher heterozygosities than PAWR. Note also that the sample "EE10D01" was assigned as females, but this is a failed sample with very few reads, and so the heterozygosity of "zero" reflects absence of data, and this is not an accurate estimation of that sample's sex.

| Heterozygosity | V1                | location.ID | location.location | location.group | location.Fst_group | location.plot_order | sex    |
|----------------|-------------------|-------------|-------------------|----------------|--------------------|---------------------|--------|
| 0.001021078    | PAWR_WIWR_AD29A07 | ED29A07     | GL                | PAWR           | PAWR               | 1                   | Male   |
| 0.001063591    | PAWR_WIWR_ED26A01 | ED26A01     | GL                | PAWR           | PAWR               | 1                   | Male   |
| 0.001127347    | PAWR_WIWR_ED27A05 | ED27A05     | GL                | PAWR           | PAWR               | 1                   | Male   |
| 0.001155423    | PAWR_WIWR_ED28A01 | ED28A01     | GL                | PAWR           | PAWR               | 1                   | Male   |
| 0.001133109    | PAWR_WIWR_ED29A02 | ED29A02     | GL                | PAWR           | PAWR               | 1                   | Male   |
| 0.002051267    | PAWR_WIWR_ED29A04 | ED29A04     | GL                | Hybrid         | Hybrid             | 1                   | Male   |
| 0.001012325    | PAWR_WIWR_ED29A05 | ED29A05     | GL                | PAWR           | PAWR               | 1                   | Male   |
| 0.001933696    | PAWR_WIWR_EE05D01 | EE05D01     | SL                | WIWR           | WIWR               | 1                   | Male   |
| 0.001646015    | PAWR_WIWR_EE05D03 | EE05D03     | SL                | WIWR           | WIWR               | 1                   | Male   |
| 0.001832412    | PAWR_WIWR_EE06D01 | EE06D01     | SL                | WIWR           | WIWR               | 1                   | Male   |
| 0.001635459    | PAWR_WIWR_EE06D02 | EE06D02     | SL                | WIWR           | WIWR               | 1                   | Male   |
| 0              | PAWR_WIWR_EE10D01 | EE10D01     | TR                | low_read_WIWR  | low_read_WIWR      | 1                   | Female |
| 0.000918246    | PAWR_WIWR_EE25D03 | EE25D03     | WH                | PAWR           | PAWR               | 1                   | Male   |
| 0.001000832    | PAWR_WIWR_EE29D02 | EE29D02     | WH                | PAWR           | PAWR               | 1                   | Male   |
| 0.001676085    | PAWR_WIWR_EF18D01 | EF18D01     | TR                | WIWR           | WIWR               | 1                   | Male   |
| 0.001027742    | PAWR_WIWR_EF18D03 | EF18D03     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001768839    | PAWR_WIWR_EF18D04 | EF18D04     | TR                | WIWR           | WIWR               | 1                   | Male   |
| 0.000972355    | PAWR_WIWR_EF22D01 | EF22D01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001045241    | PAWR_WIWR_EF23D01 | EF23D01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001170995    | PAWR_WIWR_EH01D01 | EH01D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001252033    | PAWR_WIWR_EH03D01 | EH03D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001314386    | PAWR_WIWR_EH03D02 | EH03D02     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.000198028    | PAWR_WIWR_EJ29A01 | EJ29A01     | SP                | PAWR           | PAWR               | 1                   | Female |
| 0.001371245    | PAWR_WIWR_EJ29A02 | EJ29A02     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001133101    | PAWR_WIWR_EJ29A04 | EJ29A04     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001098659    | PAWR_WIWR_EJ29A05 | EJ29A05     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.000957831    | PAWR_WIWR_EJ29A06 | EJ29A06     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001110755    | PAWR_WIWR_EK11A01 | EK11A01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.000228836    | PAWR_WIWR_EK13A02 | EK13A02     | SP                | PAWR           | PAWR               | 1                   | Female |
| 0.000415376    | PAWR_WIWR_EK13A03 | EK13A03     | SP                | PAWR           | PAWR               | 1                   | Female |
| 0.001400997    | PAWR_WIWR_FC18D01 | FC18D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001072719    | PAWR_WIWR_FD17D01 | FD17D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.00146004     | PAWR_WIWR_FD24D01 | FD24D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.00124636     | PAWR_WIWR_FD26D01 | FD26D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.00102087     | PAWR_WIWR_FD26D03 | FD26D03     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001086426    | PAWR_WIWR_FD28D01 | FD28D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.000887149    | PAWR_WIWR_FE02D03 | FE02D03     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.001541124    | PAWR_WIWR_FE13T02 | FE13T02     | KA                | PAWR           | PAWR               | 1                   | Male   |
| 0.001388576    | PAWR_WIWR_FE19T01 | FE19T01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001411893    | PAWR_WIWR_FE19T02 | FE19T02     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001033096    | PAWR_WIWR_FE20T01 | FE20T01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.000982968    | PAWR_WIWR_FE20T02 | FE20T02     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001172508    | PAWR_WIWR_FE20T03 | FE20T03     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.00119247     | PAWR_WIWR_FE22T01 | FE22T01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.000991792    | PAWR_WIWR_FE24T01 | FE24T01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001168949    | PAWR_WIWR_FE24T02 | FE24T02     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001189408    | PAWR_WIWR_FE25T01 | FE25T01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001727002    | PAWR_WIWR_FE27D01 | FE27D01     | TR                | WIWR           | WIWR               | 1                   | Male   |
| 0.001672465    | PAWR_WIWR_FF02T01 | FF02T01     | SL                | WIWR           | WIWR               | 1                   | Male   |
| 0.001094809    | PAWR_WIWR_FF07T01 | FF07T01     | HI                | PAWR           | PAWR               | 1                   | Male   |
| 0.001095629    | PAWR_WIWR_FF08T01 | FF08T01     | HI                | PAWR           | PAWR               | 1                   | Male   |
| 0.001261287    | PAWR_WIWR_FF08T02 | FF08T02     | HI                | PAWR           | PAWR               | 1                   | Male   |
| 0.001165894    | PAWR_WIWR_FF09T01 | FF09T01     | HI                | PAWR           | PAWR               | 1                   | Male   |
| 0.001613419    | PAWR_WIWR_FF11T01 | FF11T01     | SL                | WIWR           | WIWR               | 1                   | Male   |
| 0.001838493    | PAWR_WIWR_FF11T02 | FF11T02     | SL                | WIWR           | WIWR               | 1                   | Male   |
| 0.001060451    | PAWR_WIWR_FF13D01 | FF13D01     | PE                | PAWR           | PAWR               | 1                   | Male   |
| 0.001635466    | PAWR_WIWR_FF13D02 | FF13D02     | PE                | PAWR           | PAWR               | 1                   | Male   |
| 0.001162362    | PAWR_WIWR_FF15D01 | FF15D01     | CH                | PAWR           | PAWR               | 1                   | Male   |
| 0.001044182    | PAWR_WIWR_FF15D03 | FF15D03     | CH                | PAWR           | PAWR               | 1                   | Male   |
| 0.000891819    | PAWR_WIWR_FF21T01 | FF21T01     | CR                | PAWR           | PAWR               | 1                   | Male   |
| 0.001123759    | PAWR_WIWR_FF21T03 | FF21T03     | CR                | PAWR           | PAWR               | 1                   | Male   |
| 0.000965353    | PAWR_WIWR_FF22T01 | FF22T01     | CR                | PAWR           | PAWR               | 1                   | Male   |
| 0.00109898     | PAWR_WIWR_FF22T02 | FF22T02     | NE                | PAWR           | PAWR               | 1                   | Male   |
| 0.001059804    | PAWR_WIWR_FF23T01 | FF23T01     | NE                | PAWR           | PAWR               | 1                   | Male   |
| 0.001299595    | PAWR_WIWR_GD22T01 | GD22T01     | SP                | MAWR           | MAWR               | 1                   | Male   |
| 0.000983316    | PAWR_WIWR_GE05D01 | GE05D01     | SP                | PAWR           | PAWR               | 1                   | Male   |
| 0.002095563    | PAWR_WIWR_GE09D01 | GE09D01     | TR                | Hybrid         | Hybrid             | 1                   | Male   |
| 0.000916861    | PAWR_WIWR_GE11D01 | GE11D01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.000902102    | PAWR_WIWR_GE11D04 | GE11D04     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.00109209     | PAWR_WIWR_GE11D05 | GE11D05     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.000991006    | PAWR_WIWR_GE13D01 | GE13D01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.00099695     | PAWR_WIWR_GE14D01 | GE14D01     | TR                | PAWR           | PAWR               | 1                   | Male   |
| 0.000360168    | PAWR_WIWR_evl295  | evl295      | MA                | WIWR           | WIWR               | 1                   | Female |
| 0.001812327    | PAWR_WIWR_sar7443 | sar7443     | MN                | WIWR           | WIWR               | 1                   | Male   |
| 0.00101413     | PAWR_WIWR_svd2195 | svd2195     | WA                | PAWR           | PAWR               | 1                   | Male   |
| 0.000463828    | PAWR_WIWR_svd2382 | svd2382     | DN                | WIWR           | WIWR               | 1                   | Female |
| 0.000560409    | PAWR_WIWR_svd2383 | svd2383     | WN                | WIWR           | WIWR               | 1                   | Female |

