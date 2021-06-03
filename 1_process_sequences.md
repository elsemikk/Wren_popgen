# Processing sequencing data

This page details the bioinformatics pipeline used to obtain the genotype dataset from the raw GBS sequencing data. This includes demultiplexing and trimming the sequences, mapping to the reference genome, calling genotypes, and filtering the SNP dataset.  

**Input**: raw `.fastq` files with sequencing data for all individuals  
**Output**: `*.012NA`, `*.indv`, and `*.pos` files containing the genotype data for all samples. Note that these files are provided in this repository in the folder `PAWR_WIWR_012NA_files`.

# Step 1: Demultiplex the data

The starting data are raw fastq files.  
First, the data was demultiplexed using the perl script `GBS_demultiplexer_30base.pl` (the code for this script is at the bottom of this page). This requires as input a text file (`extras/barcodes_PAWR_WIWR.txt` which lists two columns: the first lists the sample names, and the second lists the barcode sequence for that sample.

Here is the code to run the demultiplexing step:  
```bash
#### 26 January 2017
# Starting the analysis of PAWR / WIWR GBS sequencing reads, from library prepared by Else Mikkelsen
# Analysis by Else and Darren
# Using new iMac in Darren’s office (3.8 GHz Intel Core i5, 64 GB memory)
# As a guide, using analysis notes for warbler genomics paper

#set up directories
mkdir Genomic_analysis
mkdir Genomic_analysis/PAWR_WIWR

#the fastq files are located in Genomic_analysis/PAWR_WIWR/
cd Genomic_analysis/PAWR_WIWR

# setup folders for clean data and extra info files:

mkdir clean_data

mkdir extras

cat > extras/barcodes_PAWR_WIWR.txt   (ctrl-d to finish)

####
ED26A01	ACGG
EE05D01	TGCT
EF18D03	CATA
EJ29A02	CGAG
FD17D01	GCTT
FE19T02	ATCA
FE27D01	GACG
FF13D01	CTGT
FF23T01	TCAA
svd2195	AGTCA
ED27A05	CATCG
EE05D03	ATCGA
EF18D04	TCGAA
EJ29A04	ACCTG
FD24D01	CTCAG
FE20T01	CGCTA
FF02T01	CCTGA
FF13D02	CGACT
GE05D01	ACGCT
svd2383	GCCAT
ED28A01	TGTGCA
EE06D01	TTGACA
EF22D01	AGCTGA
EJ29A05	TGGCAA
FD26D01	CTATCG
FE20T02	GCTGAA
FF07T01	TTCCGA
FF15D01	GACTCT
GE09D01	ATGGCG
svd2382	TCATGG
ED29A02	GTACGT
EE06D02	TAGGCT
EF23D01	GGCTAG
EJ29A06	CATGTA
FD26D03	ATTCGG
FE20T03	TGACCT
FF08T01	GCTACT
FF15D03	TCGGTA
GE11D01	CTGAGG
evl295	GCCTTA
ED29A04	GGTAGCA
EE10D01	GTGACCA
EH01D01	TTATGCA
EK11A01	ATTGGCA
FD28D01	TGGTACA
FE22T01	GACCTCA
FF08T02	TGTGCCA
FF21T01	TAGACCG
GE11D04	GGATTCA
sar7443	GATCCAA
ED29A05	AATTGCG
EE25D03	TCCAGGA
EH03D01	TCAGCAG
EK13A02	CAGTGCA
FE02D03	GTACCGA
FE24T01	TGTAACG
FF09T01	TACGATA
FF21T03	GTAAGCG
GE11D05	ATGCAAT
GD22T01	CCGGTAA
NoDNA_cntrl	AATGGACA
AD29A07	AGAATGCA
EE29D02	GAATAGCA
EH03D02	ATGAGACA
EK13A03	TGCCACCA
FE13T02	ATAGAGCA
FE24T02	ACTCGCCA
FF11T01	TAGGAACA
FF22T01	GATACGAA
GE13D01	GCACCTCA
NoBrcd_Cntrl	AGTGACAA
EF18D01	CAAGTAGA
EJ29A01	GCAAGAAT
FC18D01	ACCTACCG
FE19T01	CTACCACG
FE25T01	TAGAACGA
FF11T02	AGCAGTAA
FF22T02	GAACTGAA
GE14D01	ACTCCACG
####


perl tools/GBS_demultiplexer_30base.pl extras/barcodes_PAWR_WIWR.txt HI.4503.006.pawr_wiwr_EM1_R1.fastq HI.4503.006.pawr_wiwr_EM1_R2.fastq clean_data/PAWR_WIWR
# started 4:54pm.
# finished 6:40pm. Produced file sizes (R1 and R2 separately) ranging from about 420-1200 GB (except for one individual that produced very little: PAWR_WIWR_EE10D01)
```

# Step 2: trim the sequences

Next, we trimmed the sequences with Trimmomatic. First we need a list of the sample names to loop through.

```bash
#### 30 January 2018

## Trim data for quality with Trimmomatic

## 1. Make list of individuals

cd Genomic_analysis/PAWR_WIWR

awk '{print "PAWR_WIWR_"$1}' extras/barcodes_PAWR_WIWR.txt > extras/prefix.list.PAWR_WIWR.txt
```
This is the contents of the resulting file, `prefix.list.PAWR_WIWR.txt`:

```
PAWR_WIWR_ED26A01
PAWR_WIWR_EE05D01
PAWR_WIWR_EF18D03
PAWR_WIWR_EJ29A02
PAWR_WIWR_FD17D01
PAWR_WIWR_FE19T02
PAWR_WIWR_FE27D01
PAWR_WIWR_FF13D01
PAWR_WIWR_FF23T01
PAWR_WIWR_svd2195
PAWR_WIWR_ED27A05
PAWR_WIWR_EE05D03
PAWR_WIWR_EF18D04
PAWR_WIWR_EJ29A04
PAWR_WIWR_FD24D01
PAWR_WIWR_FE20T01
PAWR_WIWR_FF02T01
PAWR_WIWR_FF13D02
PAWR_WIWR_GE05D01
PAWR_WIWR_svd2383
PAWR_WIWR_ED28A01
PAWR_WIWR_EE06D01
PAWR_WIWR_EF22D01
PAWR_WIWR_EJ29A05
PAWR_WIWR_FD26D01
PAWR_WIWR_FE20T02
PAWR_WIWR_FF07T01
PAWR_WIWR_FF15D01
PAWR_WIWR_GE09D01
PAWR_WIWR_svd2382
PAWR_WIWR_ED29A02
PAWR_WIWR_EE06D02
PAWR_WIWR_EF23D01
PAWR_WIWR_EJ29A06
PAWR_WIWR_FD26D03
PAWR_WIWR_FE20T03
PAWR_WIWR_FF08T01
PAWR_WIWR_FF15D03
PAWR_WIWR_GE11D01
PAWR_WIWR_evl295
PAWR_WIWR_ED29A04
PAWR_WIWR_EE10D01
PAWR_WIWR_EH01D01
PAWR_WIWR_EK11A01
PAWR_WIWR_FD28D01
PAWR_WIWR_FE22T01
PAWR_WIWR_FF08T02
PAWR_WIWR_FF21T01
PAWR_WIWR_GE11D04
PAWR_WIWR_sar7443
PAWR_WIWR_ED29A05
PAWR_WIWR_EE25D03
PAWR_WIWR_EH03D01
PAWR_WIWR_EK13A02
PAWR_WIWR_FE02D03
PAWR_WIWR_FE24T01
PAWR_WIWR_FF09T01
PAWR_WIWR_FF21T03
PAWR_WIWR_GE11D05
PAWR_WIWR_GD22T01
PAWR_WIWR_NoDNA_cntrl
PAWR_WIWR_AD29A07
PAWR_WIWR_EE29D02
PAWR_WIWR_EH03D02
PAWR_WIWR_EK13A03
PAWR_WIWR_FE13T02
PAWR_WIWR_FE24T02
PAWR_WIWR_FF11T01
PAWR_WIWR_FF22T01
PAWR_WIWR_GE13D01
PAWR_WIWR_NoBrcd_Cntrl
PAWR_WIWR_EF18D01
PAWR_WIWR_EJ29A01
PAWR_WIWR_FC18D01
PAWR_WIWR_FE19T01
PAWR_WIWR_FE25T01
PAWR_WIWR_FF11T02
PAWR_WIWR_FF22T02
PAWR_WIWR_GE14D01
``` 

## Run the shell script trim.sh

Now, with the list of samples, we can run trimmomatic using a shell script. The script is written with Trimmomatic installed in the working directory in a directory called `tools`, this path would need to be modified to run in any other setup.

```bash
mkdir clean_data_trim

cat > trim.sh     #then paste text  (ctrl D to finish)
```
This is the contents of trim.sh:
```bash
#!/bin/bash
#script to trim data with trimmomatic
#created by KD Oct 30, 2014; modified by DI
#usage ./trim.sh

while read prefix

do

java -jar ./tools/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -threads 1 clean_data/"$prefix"_R1.fastq clean_data/"$prefix"_R2.fastq clean_data_trim/"$prefix"_R1.fastq clean_data_trim/"$prefix"_R1_unpaired.fastq clean_data_trim/"$prefix"_R2.fastq clean_data_trim/"$prefix"_R2_unpaired.fastq TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:30

done < extras/prefix.list.PAWR_WIWR.txt
```
Now run `trim.sh`:
```bash
# make executable:
chmod 755 trim.sh

./trim.sh

# (Had to install Java JDK before running above—used JDK 8 Update 162)

# started trim at 1:08pm; finished 2:25. Looks good.
```

# Step 3: map to the reference genome

After trimming, we mapped the sequences to the Ficedula reference genome, using bwa mem v0.7.17, picard-tools-1.97, and samtools-1.7

Here are instructions to install bwa and samtools:
```bash
# Getting ready for mapping step:

# Downloaded bwa-0.7.17
cd Genomic_analysis/PAWR_WIWR/tools/bwa-0.7.17
make  # this installs bwa

cd Genomic_analysis/PAWR_WIWR

mkdir sam
mkdir bam
mkdir log

# Install samtools:
cd tools/samtools-1.7    # and similarly for bcftools and htslib
./configure --prefix=/Applications/samtools
make
make install
```

To do the next steps, we divided the dataset into sets of 20 samples, and created seperate files listing each set of samples. This way each set could be run seperately in parallel. 
```bash
#### Split prefix list into 4 parts:

cat > extras/prefix.list.PAWR_WIWR_1to20.txt   (ctrl-d to finish)

PAWR_WIWR_ED26A01
PAWR_WIWR_EE05D01
PAWR_WIWR_EF18D03
PAWR_WIWR_EJ29A02
PAWR_WIWR_FD17D01
PAWR_WIWR_FE19T02
PAWR_WIWR_FE27D01
PAWR_WIWR_FF13D01
PAWR_WIWR_FF23T01
PAWR_WIWR_svd2195
PAWR_WIWR_ED27A05
PAWR_WIWR_EE05D03
PAWR_WIWR_EF18D04
PAWR_WIWR_EJ29A04
PAWR_WIWR_FD24D01
PAWR_WIWR_FE20T01
PAWR_WIWR_FF02T01
PAWR_WIWR_FF13D02
PAWR_WIWR_GE05D01
PAWR_WIWR_svd2383

cat > extras/prefix.list.PAWR_WIWR_21to40.txt   (ctrl-d to finish)

PAWR_WIWR_ED28A01
PAWR_WIWR_EE06D01
PAWR_WIWR_EF22D01
PAWR_WIWR_EJ29A05
PAWR_WIWR_FD26D01
PAWR_WIWR_FE20T02
PAWR_WIWR_FF07T01
PAWR_WIWR_FF15D01
PAWR_WIWR_GE09D01
PAWR_WIWR_svd2382
PAWR_WIWR_ED29A02
PAWR_WIWR_EE06D02
PAWR_WIWR_EF23D01
PAWR_WIWR_EJ29A06
PAWR_WIWR_FD26D03
PAWR_WIWR_FE20T03
PAWR_WIWR_FF08T01
PAWR_WIWR_FF15D03
PAWR_WIWR_GE11D01
PAWR_WIWR_evl295

cat > extras/prefix.list.PAWR_WIWR_41to60.txt   (ctrl-d to finish)

PAWR_WIWR_ED29A04
PAWR_WIWR_EE10D01
PAWR_WIWR_EH01D01
PAWR_WIWR_EK11A01
PAWR_WIWR_FD28D01
PAWR_WIWR_FE22T01
PAWR_WIWR_FF08T02
PAWR_WIWR_FF21T01
PAWR_WIWR_GE11D04
PAWR_WIWR_sar7443
PAWR_WIWR_ED29A05
PAWR_WIWR_EE25D03
PAWR_WIWR_EH03D01
PAWR_WIWR_EK13A02
PAWR_WIWR_FE02D03
PAWR_WIWR_FE24T01
PAWR_WIWR_FF09T01
PAWR_WIWR_FF21T03
PAWR_WIWR_GE11D05
PAWR_WIWR_GD22T01

cat > extras/prefix.list.PAWR_WIWR_61to79.txt   (ctrl-d to finish)

PAWR_WIWR_NoDNA_cntrl
PAWR_WIWR_AD29A07
PAWR_WIWR_EE29D02
PAWR_WIWR_EH03D02
PAWR_WIWR_EK13A03
PAWR_WIWR_FE13T02
PAWR_WIWR_FE24T02
PAWR_WIWR_FF11T01
PAWR_WIWR_FF22T01
PAWR_WIWR_GE13D01
PAWR_WIWR_NoBrcd_Cntrl
PAWR_WIWR_EF18D01
PAWR_WIWR_EJ29A01
PAWR_WIWR_FC18D01
PAWR_WIWR_FE19T01
PAWR_WIWR_FE25T01
PAWR_WIWR_FF11T02
PAWR_WIWR_FF22T02
PAWR_WIWR_GE14D01
```

Then, we downloaded the reference genome and indexed it before mapping the sequences to it.

```bash
# get Ficula albicollis ref genome from NCBI:
#The genome sequence GCF_000247815.1_FicAlb1.5_genomic.fna.gz is located in the folder Genomic_analysis/PAWR_WIWR/ref 

gunzip -k ref/GCF_000247815.1_FicAlb1.5_genomic.fna.gz

# index the genome using guidance from https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk  :

tools/bwa-0.7.17/bwa index -a bwtsw ref/GCF_000247815.1_FicAlb1.5_genomic.fna
# took 15 min.

/Applications/samtools/bin/samtools faidx ref/GCF_000247815.1_FicAlb1.5_genomic.fna
# took a few seconds
```

The alignment step was done using a script `align_PAWR_WIWR_toCOFL_1to20.sh`
```bash
cat > align_PAWR_WIWR_toCOFL_1to20.sh   #then paste text  (ctrl D to finish)
```
Here are the contents of `align_PAWR_WIWR_toCOFL_1to20.sh`:
```bash
#!/bin/bash
#script to align data with bwa, combine se and pe data with samtools and add RG for GATK
#created by KD Oct 30, 2014; modified by DEI 22 Nov 2014 and 15 May 2016 and 20 Nov 2016 and 26 Jan 2018 
#usage ./align_PAWR_WIWR_toCOFL.sh

# make sure these folders exist
clean_data="clean_data_trim"
sam="sam"
bam="bam"
lane="PAWR_WIWR"
runbarcode="PAWR_WIWR_Else" 
log="log"

# tell it where the executables are
bwa="tools/bwa-0.7.17/bwa"
picard="tools/picard-tools-1.97"
samtools="/Applications/samtools/bin/samtools"

while read prefix
do

## run bwa
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1.fastq $clean_data/"$prefix"_R2.fastq > $sam/"$prefix".sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1_unpaired.fastq > $sam/"$prefix".R1.unpaired.sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R2_unpaired.fastq > $sam/"$prefix".R2.unpaired.sam

## add read group headers, convert to bam, sort and index
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".sam O=$bam/"$prefix".bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R1.unpaired.sam O=$bam/"$prefix".R1.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R2.unpaired.sam O=$bam/"$prefix".R2.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE

## merge se and pe bam files with samtools and index
$samtools merge $bam/"$prefix".combo.bam $bam/"$prefix".bam $bam/"$prefix".R1.unpaired.bam $bam/"$prefix".R2.unpaired.bam
$samtools index $bam/"$prefix".combo.bam

done < extras/prefix.list.PAWR_WIWR_1to20.txt
```

Now run this script.
```bash
# make executable:
chmod 755 align_PAWR_WIWR_toCOFL_1to20.sh

./align_PAWR_WIWR_toCOFL_1to20.sh
# started 3:38pm 30Jan2018. Appears to be running well. Using 30% of CPU.
# finished 12:26am 31Jan2018
```
This was repeated for each set of 20 samples
```bash
## New terminal window:

cd Genomic_analysis/PAWR_WIWR

cat > align_PAWR_WIWR_toCOFL_21to40.sh   #then paste text  (ctrl D to finish)

--
#!/bin/bash
#script to align data with bwa, combine se and pe data with samtools and add RG for GATK
#created by KD Oct 30, 2014; modified by DEI 22 Nov 2014 and 15 May 2016 and 20 Nov 2016 and 26 Jan 2018 
#usage ./align_PAWR_WIWR_toCOFL.sh

# make sure these folders exist
clean_data="clean_data_trim"
sam="sam"
bam="bam"
lane="PAWR_WIWR"
runbarcode="PAWR_WIWR_Else" 
log="log"

# tell it where the executables are
bwa="tools/bwa-0.7.17/bwa"
picard="tools/picard-tools-1.97"
samtools="/Applications/samtools/bin/samtools"

while read prefix
do

## run bwa
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1.fastq $clean_data/"$prefix"_R2.fastq > $sam/"$prefix".sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1_unpaired.fastq > $sam/"$prefix".R1.unpaired.sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R2_unpaired.fastq > $sam/"$prefix".R2.unpaired.sam

## add read group headers, convert to bam, sort and index
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".sam O=$bam/"$prefix".bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R1.unpaired.sam O=$bam/"$prefix".R1.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R2.unpaired.sam O=$bam/"$prefix".R2.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE

## merge se and pe bam files with samtools and index
$samtools merge $bam/"$prefix".combo.bam $bam/"$prefix".bam $bam/"$prefix".R1.unpaired.bam $bam/"$prefix".R2.unpaired.bam
$samtools index $bam/"$prefix".combo.bam

done < extras/prefix.list.PAWR_WIWR_21to40.txt

—————

# make executable:
chmod 755 align_PAWR_WIWR_toCOFL_21to40.sh

./align_PAWR_WIWR_toCOFL_21to40.sh
# started 3:48pm, running well.
# With two windows running, using 50-70% of CPU. So will stop here.
# finished 12:11am 31Jan2018



#### 31 Jan 2018
# continuing with mapping step. Did half yesterday. Will do last half today:

## New terminal window:

cd Genomic_analysis/PAWR_WIWR

cat > align_PAWR_WIWR_toCOFL_41to60.sh   #then paste text  (ctrl D to finish)

--
#!/bin/bash
#script to align data with bwa, combine se and pe data with samtools and add RG for GATK
#created by KD Oct 30, 2014; modified by DEI 22 Nov 2014 and 15 May 2016 and 20 Nov 2016 and 26 Jan 2018 
#usage ./align_PAWR_WIWR_toCOFL.sh

# make sure these folders exist
clean_data="clean_data_trim"
sam="sam"
bam="bam"
lane="PAWR_WIWR"
runbarcode="PAWR_WIWR_Else" 
log="log"

# tell it where the executables are
bwa="tools/bwa-0.7.17/bwa"
picard="tools/picard-tools-1.97"
samtools="/Applications/samtools/bin/samtools"

while read prefix
do

## run bwa
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1.fastq $clean_data/"$prefix"_R2.fastq > $sam/"$prefix".sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1_unpaired.fastq > $sam/"$prefix".R1.unpaired.sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R2_unpaired.fastq > $sam/"$prefix".R2.unpaired.sam

## add read group headers, convert to bam, sort and index
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".sam O=$bam/"$prefix".bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R1.unpaired.sam O=$bam/"$prefix".R1.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R2.unpaired.sam O=$bam/"$prefix".R2.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE

## merge se and pe bam files with samtools and index
$samtools merge $bam/"$prefix".combo.bam $bam/"$prefix".bam $bam/"$prefix".R1.unpaired.bam $bam/"$prefix".R2.unpaired.bam
$samtools index $bam/"$prefix".combo.bam

done < extras/prefix.list.PAWR_WIWR_41to60.txt

—————

# make executable:
chmod 755 align_PAWR_WIWR_toCOFL_41to60.sh

./align_PAWR_WIWR_toCOFL_41to60.sh
# started 10:51am 31Jan2018, running well.
# finished 6:50pm 31Jan2018



## other terminal window:

cat > align_PAWR_WIWR_toCOFL_61to79.sh   #then paste text  (ctrl D to finish)

--
#!/bin/bash
#script to align data with bwa, combine se and pe data with samtools and add RG for GATK
#created by KD Oct 30, 2014; modified by DEI 22 Nov 2014 and 15 May 2016 and 20 Nov 2016 and 26 Jan 2018 
#usage ./align_PAWR_WIWR_toCOFL.sh

# make sure these folders exist
clean_data="clean_data_trim"
sam="sam"
bam="bam"
lane="PAWR_WIWR"
runbarcode="PAWR_WIWR_Else" 
log="log"

# tell it where the executables are
bwa="tools/bwa-0.7.17/bwa"
picard="tools/picard-tools-1.97"
samtools="/Applications/samtools/bin/samtools"

while read prefix
do

## run bwa
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1.fastq $clean_data/"$prefix"_R2.fastq > $sam/"$prefix".sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R1_unpaired.fastq > $sam/"$prefix".R1.unpaired.sam
$bwa mem -M ref/GCF_000247815.1_FicAlb1.5_genomic.fna $clean_data/"$prefix"_R2_unpaired.fastq > $sam/"$prefix".R2.unpaired.sam

## add read group headers, convert to bam, sort and index
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".sam O=$bam/"$prefix".bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R1.unpaired.sam O=$bam/"$prefix".R1.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE
java -Xmx2g -jar $picard/AddOrReplaceReadGroups.jar I=$sam/"$prefix".R2.unpaired.sam O=$bam/"$prefix".R2.unpaired.bam RGID=$lane RGPL=ILLUMINA RGLB=LIB."$prefix" RGSM="$prefix" RGPU=$runbarcode SORT_ORDER=coordinate CREATE_INDEX=TRUE

## merge se and pe bam files with samtools and index
$samtools merge $bam/"$prefix".combo.bam $bam/"$prefix".bam $bam/"$prefix".R1.unpaired.bam $bam/"$prefix".R2.unpaired.bam
$samtools index $bam/"$prefix".combo.bam

done < extras/prefix.list.PAWR_WIWR_61to79.txt

—————

# make executable:
chmod 755 align_PAWR_WIWR_toCOFL_61to79.sh

./align_PAWR_WIWR_toCOFL_61to79.sh
# started 12:53pm 31Jan2018
# finished 8:00pm 31Jan2018
```

# Step 4: call genotypes

After aligning the data, we used gatk-4.0.1.0 to call genotypes.
```bash
#### 1Feb2018

# downloaded GATK-4.0.1.1 and placed in tools folder

export PATH=$PATH:/path/to/folder/Genomic_analysis/PAWR_WIWR/tools/gatk-4.0.1.0

#### Haplotype calling:

# Will use same prefix list in the mapping step, except removing the two negative controls:

cat > extras/prefix.list.PAWR_WIWR_61to77.txt   # (ctrl-d to finish)

###
PAWR_WIWR_AD29A07
PAWR_WIWR_EE29D02
PAWR_WIWR_EH03D02
PAWR_WIWR_EK13A03
PAWR_WIWR_FE13T02
PAWR_WIWR_FE24T02
PAWR_WIWR_FF11T01
PAWR_WIWR_FF22T01
PAWR_WIWR_GE13D01
PAWR_WIWR_EF18D01
PAWR_WIWR_EJ29A01
PAWR_WIWR_FC18D01
PAWR_WIWR_FE19T01
PAWR_WIWR_FE25T01
PAWR_WIWR_FF11T02
PAWR_WIWR_FF22T02
PAWR_WIWR_GE14D01
###

mkdir gvcf

# need to prepare ".dict" file for reference (got an error about file type, but fixed by changing genome from .fna to .fa):
java -jar tools/picard-tools-1.97/CreateSequenceDictionary.jar  REFERENCE=ref/GCF_000247815.1_FicAlb1.5_genomic.fa OUTPUT=ref/GCF_000247815.1_FicAlb1.5_genomic.dict 

cat > call_snp_PAWR_WIWR_1to20.sh   #then paste text  (ctrl D to finish)
```
Here are the contents of `call_snp_PAWR_WIWR_1to20.sh`:

```bash
#!/bin/bash

while read prefix
do

## make gvcfs
gatk HaplotypeCaller -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -I bam/"$prefix".combo.bam -ERC GVCF --max-alternate-alleles 2 -O gvcf/"$prefix".gvcf.vcf

done < extras/prefix.list.PAWR_WIWR_1to20.txt
```

Then run `call_snp_PAWR_WIWR_1to20.sh`
```bash
# make executable:
chmod 755 call_snp_PAWR_WIWR_1to20.sh

./call_snp_PAWR_WIWR_1to20.sh
# started at 4:20pm, 1Feb2018
# finished 10:17am, 4Feb2018
```
This process was repeated for each batch of 20 samples.
```bash
## New window

cat > call_snp_PAWR_WIWR_21to40.sh   #then paste text  (ctrl D to finish)

-----------
#!/bin/bash
export PATH=$PATH:/Users/darrenirwin/Documents/Genomic_analysis/PAWR_WIWR/tools/gatk-4.0.1.0
while read prefix
do

## make gvcfs
gatk HaplotypeCaller -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -I bam/"$prefix".combo.bam -ERC GVCF --max-alternate-alleles 2 -O gvcf/"$prefix".gvcf.vcf

done < extras/prefix.list.PAWR_WIWR_21to40.txt
---------------------

# make executable:
chmod 755 call_snp_PAWR_WIWR_21to40.sh

bash
./call_snp_PAWR_WIWR_21to40.sh
# started 4:35pm, 1Feb2018
# Using about 55% of CPU. So will stop there.
# finished 9:15am 4Feb2018


#### 2Feb2018
# checked progress of above 2 windows. Going well, but taking about 3.5 hours for each individual, producing about 500 MB for each.
# using 52% CPU, so will try another window:

cd Genomic_analysis/PAWR_WIWR

cat > call_snp_PAWR_WIWR_41to60.sh   #then paste text  (ctrl D to finish)

-----------
#!/bin/bash
export PATH=$PATH:/path/to/folder/Genomic_analysis/PAWR_WIWR/tools/gatk-4.0.1.0
while read prefix
do

## make gvcfs
gatk HaplotypeCaller -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -I bam/"$prefix".combo.bam -ERC GVCF --max-alternate-alleles 2 -O gvcf/"$prefix".gvcf.vcf

done < extras/prefix.list.PAWR_WIWR_41to60.txt
---------------------

# make executable:
chmod 755 call_snp_PAWR_WIWR_41to60.sh

bash
./call_snp_PAWR_WIWR_41to60.sh
# started 1:22pm 2Feb2018
# 3 windows using about 77% of CPU. 
# finished 3:31am 5Feb2018


#### 4 Feb 2018
# checking progress of above. First two windows almost done, whereas third still has 5 individuals to do.
# Will now divide 61to77 into three groups, balanced among three windows according to how much left on each (so 7, 7, 3):

cd Genomic_analysis/PAWR_WIWR

cat > extras/prefix.list.PAWR_WIWR_61to67.txt   # (ctrl-d to finish)

###
PAWR_WIWR_AD29A07
PAWR_WIWR_EE29D02
PAWR_WIWR_EH03D02
PAWR_WIWR_EK13A03
PAWR_WIWR_FE13T02
PAWR_WIWR_FE24T02
PAWR_WIWR_FF11T01
###


cat > extras/prefix.list.PAWR_WIWR_68to74.txt   # (ctrl-d to finish)

###
PAWR_WIWR_FF22T01
PAWR_WIWR_GE13D01
PAWR_WIWR_EF18D01
PAWR_WIWR_EJ29A01
PAWR_WIWR_FC18D01
PAWR_WIWR_FE19T01
PAWR_WIWR_FE25T01
###


cat > extras/prefix.list.PAWR_WIWR_75to77.txt   # (ctrl-d to finish)

###
PAWR_WIWR_FF11T02
PAWR_WIWR_FF22T02
PAWR_WIWR_GE14D01
###


# one window (in window where 1-20 finished):

cat > call_snp_PAWR_WIWR_61to67.sh   #then paste text  (ctrl D to finish)

-----------
#!/bin/bash
export PATH=$PATH:/path/to/folder/Genomic_analysis/PAWR_WIWR/tools/gatk-4.0.1.0
while read prefix
do

## make gvcfs
gatk HaplotypeCaller -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -I bam/"$prefix".combo.bam -ERC GVCF --max-alternate-alleles 2 -O gvcf/"$prefix".gvcf.vcf

done < extras/prefix.list.PAWR_WIWR_61to67.txt
---------------------

# make executable:
chmod 755 call_snp_PAWR_WIWR_61to67.sh

bash
./call_snp_PAWR_WIWR_61to67.sh
# started about 12:10pm 4Feb2017
# finished 10:40am 5Feb2018

# another window (in window where 21-40 finished):

cat > call_snp_PAWR_WIWR_68to74.sh   #then paste text  (ctrl D to finish)

-----------
#!/bin/bash
export PATH=$PATH:/path/to/folder/Genomic_analysis/PAWR_WIWR/tools/gatk-4.0.1.0
while read prefix
do

## make gvcfs
gatk HaplotypeCaller -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -I bam/"$prefix".combo.bam -ERC GVCF --max-alternate-alleles 2 -O gvcf/"$prefix".gvcf.vcf

done < extras/prefix.list.PAWR_WIWR_68to74.txt
---------------------

# make executable:
chmod 755 call_snp_PAWR_WIWR_68to74.sh

bash
./call_snp_PAWR_WIWR_68to74.sh
# started 12:13pm 4Feb2018
# finished 10:49am 5Feb2018


#### 5 Feb 2018

# third window finished running. So will do the last set there:

cat > call_snp_PAWR_WIWR_75to77.sh   #then paste text  (ctrl D to finish)

-----------
#!/bin/bash
export PATH=$PATH:/path/to/folder/Genomic_analysis/PAWR_WIWR/tools/gatk-4.0.1.0
while read prefix
do

## make gvcfs
gatk HaplotypeCaller -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -I bam/"$prefix".combo.bam -ERC GVCF --max-alternate-alleles 2 -O gvcf/"$prefix".gvcf.vcf

done < extras/prefix.list.PAWR_WIWR_75to77.txt
---------------------

# make executable:
chmod 755 call_snp_PAWR_WIWR_75to77.sh

bash
./call_snp_PAWR_WIWR_75to77.sh
# started around 7:40am 5Feb2017
# finished 5:12pm 5Feb2017
```

Now we can combine these seperate genotype calls into a file with all SNPs and a file with all sites (including invariant sites, for calculating nucleotide diversity)

First, we made a file with the SNPs only.
```bash
# two approaches: 1) produce one gvcf file for SNPs only. 2) Produce one gvcf file for each chromosome, with all sites (including invariant)
# will do #1 first

# with this new version of GATK, this bit is different, involving an additional step: GenomicsDBImport.
# this has to be done for each chromosome separately.

## Make list:
cd Genomic_analysis/PAWR_WIWR

ls -o gvcf/*vcf | awk '{print "--variant "$8}' > extras/prefix.list.PAWR_WIWR.gvcf.list

cat > extras/PAWR_WIWR.sample_map   #(ctrl-d to finish)
###
PAWR_WIWR_AD29A07	gvcf/PAWR_WIWR_AD29A07.gvcf.vcf
PAWR_WIWR_ED26A01	gvcf/PAWR_WIWR_ED26A01.gvcf.vcf
PAWR_WIWR_ED27A05	gvcf/PAWR_WIWR_ED27A05.gvcf.vcf
PAWR_WIWR_ED28A01	gvcf/PAWR_WIWR_ED28A01.gvcf.vcf
PAWR_WIWR_ED29A02	gvcf/PAWR_WIWR_ED29A02.gvcf.vcf
PAWR_WIWR_ED29A04	gvcf/PAWR_WIWR_ED29A04.gvcf.vcf
PAWR_WIWR_ED29A05	gvcf/PAWR_WIWR_ED29A05.gvcf.vcf
PAWR_WIWR_EE05D01	gvcf/PAWR_WIWR_EE05D01.gvcf.vcf
PAWR_WIWR_EE05D03	gvcf/PAWR_WIWR_EE05D03.gvcf.vcf
PAWR_WIWR_EE06D01	gvcf/PAWR_WIWR_EE06D01.gvcf.vcf
PAWR_WIWR_EE06D02	gvcf/PAWR_WIWR_EE06D02.gvcf.vcf
PAWR_WIWR_EE10D01	gvcf/PAWR_WIWR_EE10D01.gvcf.vcf
PAWR_WIWR_EE25D03	gvcf/PAWR_WIWR_EE25D03.gvcf.vcf
PAWR_WIWR_EE29D02	gvcf/PAWR_WIWR_EE29D02.gvcf.vcf
PAWR_WIWR_EF18D01	gvcf/PAWR_WIWR_EF18D01.gvcf.vcf
PAWR_WIWR_EF18D03	gvcf/PAWR_WIWR_EF18D03.gvcf.vcf
PAWR_WIWR_EF18D04	gvcf/PAWR_WIWR_EF18D04.gvcf.vcf
PAWR_WIWR_EF22D01	gvcf/PAWR_WIWR_EF22D01.gvcf.vcf
PAWR_WIWR_EF23D01	gvcf/PAWR_WIWR_EF23D01.gvcf.vcf
PAWR_WIWR_EH01D01	gvcf/PAWR_WIWR_EH01D01.gvcf.vcf
PAWR_WIWR_EH03D01	gvcf/PAWR_WIWR_EH03D01.gvcf.vcf
PAWR_WIWR_EH03D02	gvcf/PAWR_WIWR_EH03D02.gvcf.vcf
PAWR_WIWR_EJ29A01	gvcf/PAWR_WIWR_EJ29A01.gvcf.vcf
PAWR_WIWR_EJ29A02	gvcf/PAWR_WIWR_EJ29A02.gvcf.vcf
PAWR_WIWR_EJ29A04	gvcf/PAWR_WIWR_EJ29A04.gvcf.vcf
PAWR_WIWR_EJ29A05	gvcf/PAWR_WIWR_EJ29A05.gvcf.vcf
PAWR_WIWR_EJ29A06	gvcf/PAWR_WIWR_EJ29A06.gvcf.vcf
PAWR_WIWR_EK11A01	gvcf/PAWR_WIWR_EK11A01.gvcf.vcf
PAWR_WIWR_EK13A02	gvcf/PAWR_WIWR_EK13A02.gvcf.vcf
PAWR_WIWR_EK13A03	gvcf/PAWR_WIWR_EK13A03.gvcf.vcf
PAWR_WIWR_FC18D01	gvcf/PAWR_WIWR_FC18D01.gvcf.vcf
PAWR_WIWR_FD17D01	gvcf/PAWR_WIWR_FD17D01.gvcf.vcf
PAWR_WIWR_FD24D01	gvcf/PAWR_WIWR_FD24D01.gvcf.vcf
PAWR_WIWR_FD26D01	gvcf/PAWR_WIWR_FD26D01.gvcf.vcf
PAWR_WIWR_FD26D03	gvcf/PAWR_WIWR_FD26D03.gvcf.vcf
PAWR_WIWR_FD28D01	gvcf/PAWR_WIWR_FD28D01.gvcf.vcf
PAWR_WIWR_FE02D03	gvcf/PAWR_WIWR_FE02D03.gvcf.vcf
PAWR_WIWR_FE13T02	gvcf/PAWR_WIWR_FE13T02.gvcf.vcf
PAWR_WIWR_FE19T01	gvcf/PAWR_WIWR_FE19T01.gvcf.vcf
PAWR_WIWR_FE19T02	gvcf/PAWR_WIWR_FE19T02.gvcf.vcf
PAWR_WIWR_FE20T01	gvcf/PAWR_WIWR_FE20T01.gvcf.vcf
PAWR_WIWR_FE20T02	gvcf/PAWR_WIWR_FE20T02.gvcf.vcf
PAWR_WIWR_FE20T03	gvcf/PAWR_WIWR_FE20T03.gvcf.vcf
PAWR_WIWR_FE22T01	gvcf/PAWR_WIWR_FE22T01.gvcf.vcf
PAWR_WIWR_FE24T01	gvcf/PAWR_WIWR_FE24T01.gvcf.vcf
PAWR_WIWR_FE24T02	gvcf/PAWR_WIWR_FE24T02.gvcf.vcf
PAWR_WIWR_FE25T01	gvcf/PAWR_WIWR_FE25T01.gvcf.vcf
PAWR_WIWR_FE27D01	gvcf/PAWR_WIWR_FE27D01.gvcf.vcf
PAWR_WIWR_FF02T01	gvcf/PAWR_WIWR_FF02T01.gvcf.vcf
PAWR_WIWR_FF07T01	gvcf/PAWR_WIWR_FF07T01.gvcf.vcf
PAWR_WIWR_FF08T01	gvcf/PAWR_WIWR_FF08T01.gvcf.vcf
PAWR_WIWR_FF08T02	gvcf/PAWR_WIWR_FF08T02.gvcf.vcf
PAWR_WIWR_FF09T01	gvcf/PAWR_WIWR_FF09T01.gvcf.vcf
PAWR_WIWR_FF11T01	gvcf/PAWR_WIWR_FF11T01.gvcf.vcf
PAWR_WIWR_FF11T02	gvcf/PAWR_WIWR_FF11T02.gvcf.vcf
PAWR_WIWR_FF13D01	gvcf/PAWR_WIWR_FF13D01.gvcf.vcf
PAWR_WIWR_FF13D02	gvcf/PAWR_WIWR_FF13D02.gvcf.vcf
PAWR_WIWR_FF15D01	gvcf/PAWR_WIWR_FF15D01.gvcf.vcf
PAWR_WIWR_FF15D03	gvcf/PAWR_WIWR_FF15D03.gvcf.vcf
PAWR_WIWR_FF21T01	gvcf/PAWR_WIWR_FF21T01.gvcf.vcf
PAWR_WIWR_FF21T03	gvcf/PAWR_WIWR_FF21T03.gvcf.vcf
PAWR_WIWR_FF22T01	gvcf/PAWR_WIWR_FF22T01.gvcf.vcf
PAWR_WIWR_FF22T02	gvcf/PAWR_WIWR_FF22T02.gvcf.vcf
PAWR_WIWR_FF23T01	gvcf/PAWR_WIWR_FF23T01.gvcf.vcf
PAWR_WIWR_GD22T01	gvcf/PAWR_WIWR_GD22T01.gvcf.vcf
PAWR_WIWR_GE05D01	gvcf/PAWR_WIWR_GE05D01.gvcf.vcf
PAWR_WIWR_GE09D01	gvcf/PAWR_WIWR_GE09D01.gvcf.vcf
PAWR_WIWR_GE11D01	gvcf/PAWR_WIWR_GE11D01.gvcf.vcf
PAWR_WIWR_GE11D04	gvcf/PAWR_WIWR_GE11D04.gvcf.vcf
PAWR_WIWR_GE11D05	gvcf/PAWR_WIWR_GE11D05.gvcf.vcf
PAWR_WIWR_GE13D01	gvcf/PAWR_WIWR_GE13D01.gvcf.vcf
PAWR_WIWR_GE14D01	gvcf/PAWR_WIWR_GE14D01.gvcf.vcf
PAWR_WIWR_evl295	gvcf/PAWR_WIWR_evl295.gvcf.vcf
PAWR_WIWR_sar7443	gvcf/PAWR_WIWR_sar7443.gvcf.vcf
PAWR_WIWR_svd2195	gvcf/PAWR_WIWR_svd2195.gvcf.vcf
PAWR_WIWR_svd2382	gvcf/PAWR_WIWR_svd2382.gvcf.vcf
PAWR_WIWR_svd2383	gvcf/PAWR_WIWR_svd2383.gvcf.vcf
###
```
First, we ran with just chromosome 28
```bash
# trying this command, using chromosome 28 (code in Ficedula genome: NC_021699.1

gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr28 -L NC_021699.1 --sample-name-map extras/PAWR_WIWR.sample_map

mkdir combined_vcfs

# now actually produce the combined vcf for chr 28

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr28 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.vcf
# took only a few minutes. Produced 49 MB file.

## will now filter:

################# Scripts used in processing
# This script is for filtering out sites with MQ lower than 20.0 
# (and leaving in sites without an MQ score, e.g. invariants):

cat > extras/vcf2minmq.pl    #then paste text  (ctrl D to finish)
```
Here is the contents of `vcf2minmq.pl`:

```bash
#!/bin/perl
use warnings;
use strict;
#This script cuts sites that have a MQ lower than the number specified. If there is no mapping quality, for example in invariant sites, then they are printed anyway.
#USAGE: cat inputfile.vcf | perl vcf2minmq.pl min_MQ > output.vcf
#Example usage for gzipped: zcat catSNPs.vcf.gz | perl /home/bin/vcf2minmq.pl 20 > catSNPs.filtered.vcf

my $minmq = $ARGV[0];

my $goodlines = 0;
my $cutlines = 0;

while (<STDIN>){
  chomp;
  my $line = $_;
  if ($. == 1){
    print "$line";
  }
  else{
    if ($line=~m/^#/){
      print "\n$line";
    }
    else{
      my @a = split(/\t/,$line);
      my $info = $a[7];
      my @infos = split(/;/,$info);
      foreach my $field (@infos){
          if ($field=~m/^MQ=/){
              $field=~s/MQ=//g;
              if ($field < $minmq){
		print STDERR"Cut $a[0]_$a[1] because MQ=$field\n";
		$cutlines++;
		goto END;
		}
          }
      }
      $goodlines++;
      print "\n$line";
    }
  }
  END:
}

print STDERR "There were $goodlines printed sites and $cutlines cut sites.\n"
```
Another script, `vcf2maxhet.pl`, is used to filter for heterozygosity:
```bash
### This next script is for filtering by heterozygosity:

cat > extras/vcf2maxhet.pl    #then paste text below  (ctrl D to finish)
```
Here are the contents of `vcf2maxhet.pl`:
```bash
#!/bin/perl
use warnings;
use strict;

#USAGE: cat inputfile.vcf | perl vcf2minhet.pl maximum_heterozygosity > output.vcf
#Example usage for gzipped: zcat catSNPs.vcf.gz | perl /home/bin/vcf2minhet.pl 0.6 > catSNPs.filtered.vcf

my $maxhet = $ARGV[0];

my $goodlines = 0;
my $cutlines = 0;

while (<STDIN>){
  chomp;
  my $line = $_;
  if ($. == 1){
    print "$line";
  }
  else{
    if ($line=~m/^#/){
      print "\n$line";
    }
    else{
      my $homo = 0;
      my $het = 0;
      my @a = split(/\t/,$line);
      foreach my $i(9..$#a){
        my $info = $a[$i];
        my @infos = split(/:/,$info);
        if ($infos[0] eq "./."){
          next;
        }
        my @alleles = split(/\//, $infos[0]);
        if ($alleles[0] eq $alleles[1]){
          $homo++;
        }else{
          $het++;
        }
      }
      my $total = $homo + $het;
      if ($total == 0){
        #print STDERR "Cut $a[0]_$a[1] because no data\n";
        next;
        $cutlines++;
      }
      my $hetperc = $het / $total;
      if ($hetperc < $maxhet){
        print "\n$line";
        $goodlines++;
      }
      else{
        print STDERR "Cut $a[0]_$a[1] because Het observed = $hetperc\n";
        $cutlines++;
      }
    }
  }
}

print STDERR "There were $goodlines printed sites and $cutlines cut sites.\n"
```
Now we filtered the chromosome 28 file using vcftools and the custom scripts above.
```bash
# Remove indels and SNPs with more than two alleles:
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.vcf --remove-indels --max-alleles 2 --recode --recode-INFO-all --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.vcf  
# After filtering, kept 77 out of 77 Individuals
# Outputting VCF file...
# After filtering, kept 15708 out of a possible 20951 Sites
# Run Time = 3.00 seconds

# Remove sites with more than 60% missing genotypes (note the parameter below is intuition times minus one):
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.vcf.recode.vcf --max-missing 0.4 --recode --recode-INFO-all --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.vcf
# After filtering, kept 77 out of 77 Individuals
# Outputting VCF file...
# After filtering, kept 7214 out of a possible 15708 Sites
# Run Time = 2.00 seconds

# Filter out sites with MQ lower than 20.0 (Greg's script):
cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.vcf.recode.vcf | perl extras/vcf2minmq.pl 20.0 > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.MQ20.vcf
# There were 7022 printed sites and 192 cut sites.

# Filter out sites with heterozygosity above 60% (Greg's script):
cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.MQ20.vcf | perl extras/vcf2maxhet.pl 0.6 > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.MQ20.lowHet.vcf
# There were 7010 printed sites and 12 cut sites.

# Convert to tab file in 012NA format (run two commands below):
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.MQ20.lowHet.vcf --012 --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.MQ20.lowHet.tab

cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.MQ20.lowHet.tab.012 | sed 's/-1/NA/g' > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr28.max2allele_noindel.maxmiss60.MQ20.lowHet.tab.012NA
# file is 1.2 MB in size
```
This produces the 012NA, indv, and pos files for chromosome 28 which are used for subsequent analyses.

After testing chromosome 28, we ran the next 4 chromosomes using a bash script:
```bash
cat > make_combined_vcf_for_chromosomes1to3.sh   #then paste text  (ctrl D to finish)
```
Here are the contents of `make_combined_vcf_for_chromosomes1to3.sh`:
```bash
#chr1:
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr1 -L NC_021671.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr1 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr1.vcf

#chr1A
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr1A -L NC_021672.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr1A -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr1A.vcf

#chr2
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr2 -L NC_021673.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr2 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr2.vcf

#chr3
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr3 -L NC_021674.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr3 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr3.vcf
```
Run this script:
```bash
chmod 755 make_combined_vcf_for_chromosomes1to3.sh

./make_combined_vcf_for_chromosomes1to3.sh
# started 9:35pm 5Feb2018
# finished 10:49pm same day
```
Next, we ran processed chromosomes 4-8 in a bash script:
```bash
cat > make_combined_vcf_for_chromosomes4to8.sh   #then paste text  (ctrl D to finish)
```
Here are the contents of `make_combined_vcf_for_chromosomes4to8.sh`:
```bash
#chr4
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr4 -L NC_021675.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr4 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr4.vcf

#chr4A
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr4A -L NC_021676.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr4A -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr4A.vcf

#chr5
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr5 -L NC_021677.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr5 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr5.vcf

#chr6
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr6 -L NC_021678.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr6 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr6.vcf

#chr7
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr7 -L NC_021679.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr7 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr7.vcf

#chr8
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr8 -L NC_021680.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr8 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr8.vcf
```
Then, run the script.
```bash
chmod 755 make_combined_vcf_for_chromosomes4to8.sh

./make_combined_vcf_for_chromosomes4to8.sh
# started about 9:40pm 5Feb2018
# finished 10:40pm same day
```
Next, chromosomes 9 and above were processed in serial using a bash script.
```bash
#### 6 Feb 2018
# continuing with producing combined vcf files for each chromosome

cat > make_combined_vcf_for_chromosomes9toAll.sh   #then paste text  (ctrl D to finish)
```
Here are the contents of `make_combined_vcf_for_chromosomes9toAll.sh`:
```bash
#chr9
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr9 -L NC_021681.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr9 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr9.vcf

#chr10
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr10 -L NC_021682.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr10 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr10.vcf

#chr11
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr11 -L NC_021683.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr11 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr11.vcf

#chr12
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr12 -L NC_021684.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr12 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr12.vcf

#chr13
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr13 -L NC_021685.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr13 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr13.vcf

#chr14
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr14 -L NC_021686.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr14 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr14.vcf

#chr15
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr15 -L NC_021687.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr15 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr15.vcf

#chr17
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr17 -L NC_021688.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr17 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr17.vcf

#chr18
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr18 -L NC_021689.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr18 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr18.vcf

#chr19
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr19 -L NC_021690.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr19 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr19.vcf

#chr20
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr20 -L NC_021691.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr20 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr20.vcf

#chr21
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr21 -L NC_021692.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr21 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr21.vcf

#chr22
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr22 -L NC_021693.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr22 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr22.vcf

#chr23
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr23 -L NC_021694.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr23 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr23.vcf

#chr24
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr24 -L NC_021695.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr24 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr24.vcf

#chr25
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr25 -L NC_021696.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr25 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr25.vcf

#chr26
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr26 -L NC_021697.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr26 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr26.vcf

#chr27
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chr27 -L NC_021698.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chr27 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr27.vcf

#chrZ
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chrZ -L NC_021700.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chrZ -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chrZ.vcf

#chrLG34
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chrLG34 -L NC_021701.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chrLG34 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chrLG34.vcf

#chrLG35
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chrLG35 -L NC_021702.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chrLG35 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chrLG35.vcf

#chrLGE22
gatk GenomicsDBImport --genomicsdb-workspace-path PAWR_WIWR_combined_genomicsDB_chrLGE22 -L NC_021703.1 --sample-name-map extras/PAWR_WIWR.sample_map

gatk GenotypeGVCFs -R ref/GCF_000247815.1_FicAlb1.5_genomic.fna -V gendb://PAWR_WIWR_combined_genomicsDB_chrLGE22 -O combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chrLGE22.vcf
```
Then, run the script:
```bash
chmod 755 make_combined_vcf_for_chromosomes9toAll.sh

./make_combined_vcf_for_chromosomes9toAll.sh
# started 9:31am 6Feb2018
# finished 10:59am same day
```

Now that all chromosomes were processed into single files, we filtered them using vcftools (as was done previously for chr28).

To do this, we looped through each chromosome using a text file listing the names of the chromosomes:
```bash
cat > extras/chromosomes_to_process_all.txt    #then paste text  (ctrl D to finish)
```
Here are the contents of `chromosomes_to_process_all.txt`:
```
1
1A
2
3
4
4A
5
6
7
8
9
10
11
12
13
14
15
17
18
19
20
21
22
23
24
25
26
27
Z
LG34
LG35
LGE22
```

We used a bash script to process the chromosomes.
```bash
cat > filter_and_prep_PAWR_WIWR_by_chr.sh   #then paste text  (ctrl D to finish)
```
Here are the contents of `filter_and_prep_PAWR_WIWR_by_chr.sh`:
```bash
#!/bin/bash

while read chr
do

# Remove indels and SNPs with more than two alleles:
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".vcf --remove-indels --max-alleles 2 --recode --recode-INFO-all --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.vcf  

# Remove sites with more than 60% missing genotypes (note the parameter below is intuition times minus one):
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.vcf.recode.vcf --max-missing 0.4 --recode --recode-INFO-all --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.vcf

# Filter out sites with MQ lower than 20.0 (Greg's script):
cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.vcf.recode.vcf | perl extras/vcf2minmq.pl 20.0 > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.MQ20.vcf

# Filter out sites with heterozygosity above 60% (Greg's script):
cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.MQ20.vcf | perl extras/vcf2maxhet.pl 0.6 > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.MQ20.lowHet.vcf

# Convert to tab file in 012NA format (run two commands below):
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.MQ20.lowHet.vcf --012 --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.MQ20.lowHet.tab

cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.MQ20.lowHet.tab.012 | sed 's/-1/NA/g' > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss60.MQ20.lowHet.tab.012NA

done < extras/chromosomes_to_process_all.txt
```
Run this script:
```bash
chmod 755 filter_and_prep_PAWR_WIWR_by_chr.sh

./filter_and_prep_PAWR_WIWR_by_chr.sh
# started about 5:18pm 5Feb2018
# finished 5:32pm same day
```
This produced the `*.012NA, *.indv, and *.pos` files for further analysis.

We also produced another dataset with stricter missing data requirements (max 30% missing) using the same workflow as above. The analyses in the paper were performed using this stricter missing data cutoff.

```bash
## Will make another set with SNPs removed if missing more than 30% of genotypes (because PCA not so good with lots of missing data):

cat > extras/chromosomes_to_process_all.txt    #then paste text  (ctrl D to finish)

#####
1
1A
2
3
4
4A
5
6
7
8
9
10
11
12
13
14
15
17
18
19
20
21
22
23
24
25
26
27
28
Z
LG34
LG35
LGE22
#####

cat > filter_maxmiss30_and_prep_PAWR_WIWR_by_chr.sh   #then paste text  (ctrl D to finish)

#####
#!/bin/bash

while read chr
do

# Remove sites with more than 30% missing genotypes (note the parameter below is intuition times minus one):
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.vcf.recode.vcf --max-missing 0.7 --recode --recode-INFO-all --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.vcf

# Filter out sites with MQ lower than 20.0 (Greg's script):
cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.vcf.recode.vcf | perl extras/vcf2minmq.pl 20.0 > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.MQ20.vcf

# Filter out sites with heterozygosity above 60% (Greg's script):
cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.MQ20.vcf | perl extras/vcf2maxhet.pl 0.6 > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.MQ20.lowHet.vcf

# Convert to tab file in 012NA format (run two commands below):
tools/vcftools_0.1.12b/bin/vcftools --vcf combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.MQ20.lowHet.vcf --012 --out combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.MQ20.lowHet.tab

cat combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.MQ20.lowHet.tab.012 | sed 's/-1/NA/g' > combined_vcfs/PAWR_WIWR.genotypes.SNPs_only.chr"$chr".max2allele_noindel.maxmiss30.MQ20.lowHet.tab.012NA

done < extras/chromosomes_to_process_all.txt
#####

chmod 755 filter_maxmiss30_and_prep_PAWR_WIWR_by_chr.sh

./filter_maxmiss30_and_prep_PAWR_WIWR_by_chr.sh
```
This produced the `*.012NA, *.indv, and *.pos` files for analysis.
