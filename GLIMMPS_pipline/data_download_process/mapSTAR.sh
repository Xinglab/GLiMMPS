#!/bin/bash
##
## Align the reads from pairwise fastq files to hg19 genmoe using STAR and
## get the uniquely mapped reads
## Author: Yungang Xu (yungang.xu@hotmail.com)
##----------Arguments----------
## $1 the individual id, which is used to find the pairwise fastq files for mapping
##----------Usage----------
## bash mapSTAR.sh <Individual\_ID>

## set the input, output, and STAR index path
inPrefix='/mnt/isilon/xing_lab/aspera/Yungang/Geuvadis/CEU/input/' ## the path prefix of the input fastq files
outPrefix='/mnt/isilon/xing_lab/xuy5/GLIMMPS/Geuvadis/CEU/' ## the path prefix of the output files
hg19_starIndex_75='/mnt/isilon/xing_lab/xuy5/Genome/hg19/starIndex' ## the parth of the STAR index for hg19 genome

## do mapping
echo STAR mapping
STAR --runThreadN 50 --genomeDir $hg19_starIndex_75 --readFilesCommand gunzip -c --outFileNamePrefix $outPrefix'STARout/'$1. --readFilesIn $inPrefix$1'_1.fq.gz' $inPrefix$1'_2.fq.gz' --outFilterMultimapNmax 20 --alignEndsType EndToEnd --outFilterMismatchNmax 6 --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignIntronMax 300000 --outSJfilterOverhangMin 8 8 8 8 SortedByCoordinate --outFilterMultimapScoreRange 0
rm -r $outPrefix'STARout/'$1'._STARtmp'
mv $outPrefix'STARout/'$1.Aligned.out.sam $outPrefix'sam/'$1.sam

## get uniquely mapped reads
echo Get uniquely mapped reads
mkdir -p $outPrefix/sam_unique/$1
## output the SAM header
samtools view -S -H $outPrefix/sam/$1.sam > $outPrefix/sam_unique/$1.unique.sam
samtools view -S $outPrefix/sam/$1.sam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^NH:i:"){ mp = substr($i,6,length($i)-5); if (mp == 1){print $0}} }}' >> $outPrefix/sam_unique/$1/unique.sam
echo $1	Done >> $1.mapSTAR.log
