#!/bin/bash
##
## The purpose of this code is to change the file format from .vcf to .tped and .tfam for downstream calculation
## install vcftools first
## Author: Yungang Xu (yungang.xu@hotmail.com)
##----------Arguments----------
## $1: chromosome number (1 to 22)
## $2: population id (CEU)
## $3: sample list file (CEU.samples.txt)
## Command parameter:
## --gzvcf: read compressed (gzipped) VCF files directly
## --plink-tped: outputs the genotype data in the PLINK transposed format with suffixes ".tped" and ".tfam"
## --keep: Provide files containing a list of individuals to include in subsequent analysis.
##--------Usage----------
## bash vcf2plink.sh 1 CEU CEU.samples.txt

vcf_path='path/to/folder/of/all/vcf_files/' ## change this accordingly
samples='path/to/Geuvadis/'$2'/'$3  ## change this accordingly
mkdir -p SNPs_plink
vcftools --gzvcf /mnt/isilon/xing_lab/aspera/PAIRADISE/YG_test/allVCFs/ALL.chr$1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --plink-tped --out SNPs_plink/$2.plink.$1 --keep $samples

