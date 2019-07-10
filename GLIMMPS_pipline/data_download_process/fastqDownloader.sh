#!/bin/bash
## Author: Yungang Xu (yungang.xu@hotmail.com)
## Batch download Geuvadis fastq files
##-----------Arguments----------
## $1 the metadata (sample list) downloaded from http://www.internationalgenome.org/data-portal/data-collection/geuvadis
## remove the header line before use it
## $2 a number indicating the parallel threads for multiple downloading simultaneously. Based on my experience,
##    the above server could afford upto 70 threads, otherwise will get some download failed because of the connection traffic.
##----------USAGE----------
## [qsub -cwd -l h_vmem=120G,m_mem_free=60G -M xuy5@email.chop.edu -m bea] fastqDownloader.sh igsr_Geuvadis.tsv.tsv 20

download () {
	url=$1 # the ftp url
	indi=$6 # the individual ID, same with vcf file
	pop=$7 # population ID
	mkdir -p $pop # group data based on their population ID
	fqid=$(basename $url) # get the fqid, i.e. 1.fastq.gz or 2.fastq.gz
	wget -c $url -O ${pop}/${indi}_${fqid}
	echo $@ Done >> download.${pop}.log
}
a=0
while read -r line
do
if [ $a -lt $2 ]
then
	download $line &
	a=`expr $a + 1`
else
	download $line
	a=0
fi
done < $1
