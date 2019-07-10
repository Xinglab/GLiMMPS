#!/usr/bin/perl -w

##　For all individuals in the population, obtain the unique reads from the tophat mapping result file: accepted_hits.bam.

open ( IN, "<CheungCEU_IDs.txt");  ## The file list all Indivividual IDs in the population.
my @ceus = <IN> ;

open (OUT, ">submit_uniquereads");
foreach my $ceuid (@ceus) {
chomp($ceuid) ;
print $ceuid,"\n" ;

#my $awkstr = q (awk  -F"\t" '($2~"P"&&($6=="50M"||($6~"N"&&$6!~"D"&&$6!~"I"))) || NF<7') ; 
my $awkstr = q (awk  -F"\t" '($5=="255"  &&($6=="50M"||($6~"N"&&$6!~"D"&&$6!~"I"))) || NF<7') ; 
print OUT "
samtools view -hX -q 255 -o $ceuid/temp.sam $ceuid/accepted_hits.bam
$awkstr $ceuid/temp.sam > $ceuid/unique.sam

" ;

}
close(OUT) ;
#
system("chmod +x submit_uniquereads");
system("./submit_uniquereads");
