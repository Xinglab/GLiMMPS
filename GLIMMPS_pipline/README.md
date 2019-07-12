# GLiMMPS_pipeline
## About
This repository contains the source codes of the pipeline for sQTL mapping and illustrates the usage by performing an example analysis on the data of Gevadis CEU population.
## Table of Contents
- [1. Prerequisites](#1)
  - [1.1 STAR](#1.1)
  - [1.2 samtools](#1.2)
  - [1.3 vcftools](#1.3)
  - [1.4 plink](#1.4)
- [2. Data preparation](#2)
  - [2.1 RNA-seq data downlaod](#2.1)
  - [2.2 QC and Alignment](#2.2)
  - [2.3 Genotype data](#2.3)
- [3. Run GLiMMPs](#3)
  - [3.1 Get AS event and junction read counts](#3.1)
  - [3.2 Run statistic models for sQTLs](#3.2)
- [Contacts](#4)
- [Citation](#5)

## <a name="1"></a>1. Prerequisites
The full pipeline starts with data download and ends with sQTL mapping using GLiMMPs, which requires all tools for fastq files download, quality control (QC), alignment, vcf manipulation, and running environments for GLiMMPs (Python, Perl, and R).
### <a name="1.1"></a>1.1 Download and install RNA-seq reads alignment tool `STAR` (version 2.6.1)
Refer to [STAR GitHub](https://github.com/alexdobin/STAR) page for detailed instructions.
### <a name="1.2"></a>1.2 Download and install BAM file tool `samtools`(version 1.9)
Refer to [samtools GitHub](https://github.com/samtools/samtools) page for detailed instructions.
### <a name="1.3"></a>1.3 Download and install `vcftools` (version 0.1.13)
Refer to [vcftools](https://vcftools.github.io/examples.html) page for detailed instructions.
### <a name="1.4"></a>1.4 Downlaod and install `plink` (version 1.07)
Refer to [plink](http://zzz.bwh.harvard.edu/plink/download.shtml) page for detailed instructions.

## <a name="2"></a>2. Data download and processing (Geuvadis as example)
All shell scripts (.sh files) mentioned hereafter are located in the folder *data_download_process/*.
### <a name="2.1"></a>2.1 RNA-seq data download and processing
**(1) Download metadata** (sample list): Go to Geuvadis [Dataport](http://www.internationalgenome.org/data-portal/data-collection/geuvadis), select the *Data types* as `Sequence` and the *Analysis groups* as `mRNA`, and then click `Download the list` button to get the metadata (*igsr_Geuvadis.tsv.tsv*) of all five Geuvadis populations (CEU, FIN, GBR, TSI, and YRI).

**(2) Batch download**: use `fastqDownloader.sh` to downlaod the samples you need. Take CEU population as example, first take out those lines for CEU from *igsr_Geuvadis.tsv.tsv*, you can run `grep "CEU" igsr_Geuvadis.tsv.tsv >> CEU_igsr_Geuvadis.tsv.tsv`; 

then run `bash fastqDownloader.sh CEU_igsr_Geuvadis.tsv.tsv 50` to download the fastq files for samples included in *CEU_igsr_Geuvadis.tsv.tsv* with 50 threads parallelly.

fastq files will be downloaded and grouped into different folders named as the population ID, like CEU etc. All fastq files will also be renamed as <Idividual_ID>\_<Sample_ID>\_1.fastq.gz, like NA06984_ERR188325_1.fastq.gz.

### <a name="2.2"></a>2.2 QC and alignment
**(1) QC**: we ran trim_galore on the pairwise fastq files (set --paired) to make sure all reads have the same length of 75bp; and then run fastqc for QC filtering.

**(2) Merge replicates for single individual**: after trim and QC, many of the idviduals will have more than one replicates, we merged the fastq files of replicates, so that each individual will have only one pair fastq files, named as <Individual_ID>\_1.fq.gz and <Individual_ID>\_2.fq.gz. In total, for CEU population, there will be 184 fastq files for 92 individuals.

**(3) Alignment using STAR**: run `mapSTAR.sh <Individual_ID>` to map the reads to hg19 genome and get the uniquely mapped reads. The unique.sam files will be in *sam_unique/* with the subfolder named with the <Individual_ID>.
### <a name="2.3"></a>2.3 Genotypes (vcf files) download and processin
**(1) Download vcf files**: run `vcfDownloader.sh` to download the vcf files from 1000 Genome data portal.

**(2) convert the vcf files to .tped and .tmap**: The downstream sQTL analysis will use genotypes from `plink` formatted files with each chromosome in one file. Run `vcf2plink.sh` to get the formatted files from vcf files. Output files will be in **SNPs_plink/**.

## <a name="3"></a>3. Run GLiMMPs pipeline
### <a name="3.2"></a>3.1 Obtain the junction read counts for all possible alternative splicing events in the population
**(1) Get all possible alternative splicing (AS) events in the population.**

Customize the config.GLiMMPS.txt for the population, which contains main parapeters for all data processing.

  ```./GLiMMPscode/pythonperlsrcs/batch_allASevents.pl /GLiMMPscode/config.GLiMMPS.txt```

The above code will produce a job file named **submit_ASevents**, which has the command line to run `/GLiMMPScode/pythonperlsrcs/processGTF.SAMs.py` on the given population data. You can run it directly or submit this file to a HPC.

The resulting AS events will be outputted in **/ASEvents/**.

**(2) Get junction read count for all individuals**: 

run the following code to get the job file and configure files:
```./GLiMMPscode/pythonperlsrcs/batch_getASreads.pl /GLiMMPscode/config.GLiMMPS.txt```
The above code will produce a folder named **tempconfig/**, which has the configurations to run embeded rMATs code for each individual of the given population data. 

The above code will also produce a job file named **submit_ASreadcounts**, which has the command lines to run `/GLiMMPScode/pythonperlsrcs/rMATS.processsUnique.singlesam.py` for each individual. You can run it one by one directly or submit this file to a HPC to run parallelly.

After all above jobs were done, the resulted junction read counts for all 4 types of AS events per individual will be output in **myOutput/**. Then run the following code to sumerize the read counts across individuals:
```./GLiMMPscode/pythonperlsrcs/summarizeallexoninc.pl /GLiMMPscode/config.GLiMMPS.txt```
The above code will put the summerized read counts into **Exon_Inc_Simple/AScounts/**. 

**(3) Filter the AS events**: 
Finally run the following code to filter out the AS events with no or little change in exon inclusion level (|Δψ|<0.1) or few total junction read counts (median n <5) in the population.
```
cd Exon_Inc_Simple
./GLiMMPScode/Rscripts/summarystat_exonmin5.R
```
The above code will produce files in **Exon_Inc_Simple/alltype/**, which will be used as input for statictical models.

### <a name="3.2"></a>3.2 Run statistical models for sQTLs
After above step, we have the AS exon information, the inclusion junction and total junction read counts matrix for all individuals (in **Exon_Inc_Simple/alltype/**), as well ase the plink-format genotype files (in **SNPs_plink/**).

The sQTL analysis can be performed for each chromosome using script in **/GLiMMPScode/run.GLiMMPS.parallel.sh**.

**Parameters/arguments for this script:**

  - **chr_id**  \<INT\>, the chromosome to be analyzed, number only without 'chr', i.e. {1..22}
  - **population_id** \<STR\>, the population ID, i.e. {CEU, GBR, ...}, which will be used to locate the genotype files
  - **path_to_AS**  \<STR\>, the path to **alltype/** generated by step 3.2.
  - **path_to_genotype**  \<STR\>, path to **SNPs_plink/** generated by step 2.3.
  - **thread_number** \<INT\>, the mumber of threads to used for parallel running. Defalt = 100 (required ~50GB memory for CEU).
  
**Usage:** run statistic models on chr1 for CEU population with 100 threads.

  ```./runGLiMMPS.parallel.sh 1 CEU alltype SNPS_plink 100```
  
**Outputs:**
  
 Two folders: 
  - **parallel.tmpasso/** includes all associations, one file per AS event; one subfolder per chromosome; and one PDF per AS event that has at least one significant sQTL. These are the main results of sQTLs.
  - **parallel.tmpgeno/** includes all tempory files used by the statistical analysis. For debugging only.

## <a name="4"></a>Contacts
Yi Xing (Yi.Xing@pennmedicine.upenn.edu)

Yungang Xu (yungang.xu@hotmail.com)

## <a name="5"></a>Citation
\[1\] Zhao, K., Lu, Z. X., Park, J. W., Zhou, Q., & Xing, Y. (2013). GLiMMPS: robust statistical model for regulatory variation of alternative splicing using RNA-seq data. Genome biology, 14(7), R74.
