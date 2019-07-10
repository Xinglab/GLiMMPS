#!/bin/bash
## batch download the vcf files from 1000 Genome ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
## one vcf.gz file for each chromosome
##
## Author: Yungang Xu (yungang.xu@hotmail.com)
##--------------Arguments-----------------
## NONE
##--------------Usage-----------------
## bash vcfDownloader.sh

for i in {1..22}; 
do  
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz; 
done
