#!/bin/bash
## merge the replicates for each individual
## Authour: Yungang Xu (yungang.xu@hotmail.com)
##----------Arguments----------
## $1 individual list file, one ID per line.
##----------Usgae----------
## bash repMerge.sh CEU.samples.txt
## Output: merged/<Individual_ID>_1.fq.gz and merged/<Individual_ID>_2.fq.gz for each individual with ID <Individual_ID>.

doMerge(){
  fq1=`ls $1*_1*`
  fq2=`ls $1*_2*`
  cat $fq1 > merged/$1'_1.fq.gz'
  cat $fq2 > merged/$1'_2.fq.gz'
  echo $1	$fq1 >> merge.log.txt
  echo $1	$fq2 >> merge.log.txt
}

while read -r line
do
  doMerge $line
done < $1
