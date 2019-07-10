#!/bin/bash
## Author: Yungang Xu (yungang.xu@hotmail.com)
##-------------Arguments-----------
## $1 chromosome id (without chr) 1,2,...,22
## $2 population ID, CEU, GBR....
## $3 path/to/alltype/, the path where the alltype/ folder located
## $4 path/to/genotype/, the path where the SNPs_plink/ folder located
## $5 [optional] the number threads for parallel running; defualt = 100. A bigger number may require more memory (100 threads need ~50GB).
##-------------Usage-----------
## bash runGLiMMPS.paallel.sh 1 CEU [100]
##-------------Output----------
## two folders:
## parallel.tmpasso/ includes all associations, one file per AS event; one subfolder per chromosome; and one PDF per AS event that has at 
##                   least one significant sQTL. This the main results of sQTLs
## parallel.tmpgeno/ includes all tempory files used by the statistical analysis. Usually for debugging only.

## You may need to change the following paths so that the files could be fund.
exonInfo=$3/bychrs/exonsinfor.plink.5reads.chr$1.txt
allJunction=$3/bychrs/plink.5reads.allreads.chr$1.txt
InclJunction=$3/bychrs/plink.5reads.IJ.chr$1.txt
genotype=$4/SNPs_plink/$2.plink.$1

count=$(cat $exonInfo | wc -l ) # get the number of AS events
for (( i=1; i<$count; i++ )) # run GLiMMPS on all AS events
do
echo Processing AS $i
echo 
Rscript /mnt/isilon/xing_lab/xuy5/GLIMMPS/Geuvadis/GLiMMPScode/Rscripts/sQTLregress.oneexon.parallel.R $exonInfo $allJunction $InclJunction $genotype $i $5
done
