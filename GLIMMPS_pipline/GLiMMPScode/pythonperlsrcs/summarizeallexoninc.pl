#!/usr/bin/perl -w

####################
# summarize AS events for all individuals for all 4 types of AS events.


######################################################################
# read in the configuration file 
if ( ($#ARGV+1) != 1 ) {
	print "Usage: batch_allASevents.pl configfile \n";
	exit;
}
my $configfile =$ARGV[0];



#######################
# read in the configuration file 
#

# usage: readConfig("filename")
# returns: num. of parameters set
my @config_keys = (); 

sub readConfig
{
    my $filename=shift;
    my $i=0;
    
    if(open(CONF, "<$filename"))
    {
        while(<CONF>)
        {
	 # print $_ ;
            if(/^\s*(\S+)\s*=\s*([^(#*|\s*)]+)/)
            {
                $$1=$2;
                ++$i;
		push @config_keys, $1 ;
            }
        }
        return($i);
    }
    else
    {
        return(0);
    }
}


###########
#  
my $numparams = readConfig($configfile);
print "Read $numparams parameters from $configfile:\n\n" ;

foreach $Config_key (@config_keys) {
  
  print "$Config_key = $$Config_key\n";
  
}
######################################################################



open ( IN0, "<$IDFILE");
my @ids = <IN0> ;
close(IN0) ;

my $id1 = $ids[0];
chomp($id1) ;

#make directory unless it already exists
mkdir "Exon_Inc_Simple", 0777 unless -d "Exon_Inc_Simple";
mkdir "Exon_Inc_Simple/AScounts", 0777 unless -d "Exon_Inc_Simple/AScounts";


my @types = ("MXE","SE","A5SS","A3SS"); #,"RI");
foreach my $type (@types) # ="SE" ;
  {
my $file1 = "myOutput/JC.${id1}.$type.MATS.input.txt"; 

# Add splicing type to last column of each file


$cmmd = qq(perl -nle 'print "${type}_\$_"'  ASEvents/$PROJECTTITLE.$type.txt  >  ASEvents/${PROJECTTITLE}2.$type.txt);
print $cmmd,"\n";
system($cmmd) ;

`cut -f1 ASEvents/${PROJECTTITLE}2.$type.txt  > Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.IJ.txt`;
`cut -f1 ASEvents/${PROJECTTITLE}2.$type.txt  > Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.SJ.txt`;
`cut -f4,5 $file1 |paste  ASEvents/${PROJECTTITLE}2.$type.txt  - > Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.exoninfor.txt`;



#`cut -f1-7 $file1 > Exon_Inc_Simple/AScounts/$PROJECTTILE_PSI.txt`;
#`grep Exon_ID_list Exon_Inc_Simple/AScounts/$PROJECTTILE_UJ.txt >Exon_Inc_Simple/AScounts/title.txt`;



my $idstr ="" ;
foreach my $id (@ids) {
chomp($id) ;
#print $id,"\n" ;
$idstr .="\t$id" ;

my $file = "myOutput/JC.${id}.$type.MATS.input.txt"; 

`cut -f2 $file |sed -e 's/IJC_//' |paste Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.IJ.txt - > Exon_Inc_Simple/AScounts/tmp.txt`;
`mv Exon_Inc_Simple/AScounts/tmp.txt Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.IJ.txt `;
wait ;

`cut -f3 $file |sed -e 's/SJC_//' |paste Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.SJ.txt - > Exon_Inc_Simple/AScounts/tmp.txt`;
`mv Exon_Inc_Simple/AScounts/tmp.txt Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.SJ.txt `;
wait ;

}
if ($type eq "MXE") {
`cat Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.exoninfor.txt > Exon_Inc_Simple/AScounts/$PROJECTTITLE.alltype.exoninfor.txt`;
`cat Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.IJ.txt > Exon_Inc_Simple/AScounts/$PROJECTTITLE.alltype.IJ.txt`;
`cat Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.SJ.txt > Exon_Inc_Simple/AScounts/$PROJECTTITLE.alltype.SJ.txt`;
}

else {
`cat Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.exoninfor.txt |grep -v ${type}_ID >> Exon_Inc_Simple/AScounts/$PROJECTTITLE.alltype.exoninfor.txt`;
`cat Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.IJ.txt |grep -v ${type}_ID >> Exon_Inc_Simple/AScounts/$PROJECTTITLE.alltype.IJ.txt`;
`cat Exon_Inc_Simple/AScounts/$PROJECTTITLE.$type.SJ.txt |grep -v ${type}_ID >> Exon_Inc_Simple/AScounts/$PROJECTTITLE.alltype.SJ.txt`;
}


}
#`echo $idstr >> Exon_Inc_Simple/title.txt`;





#### Notice that MXE type have 2 extra columns for the 2nd exon. After this script, need to manually add 2 extra empty columns for the other types in file: Exon_Inc_Simple/AScounts/$PROJECTTITLE.alltype.exoninfor.txt

my $OUTDIR = "Exon_Inc_Simple/AScounts/";

 $cmmd = qq(
grep -v ^MXE $OUTDIR/$PROJECTTITLE.alltype.exoninfor.txt |cut -f1-7 >$OUTDIR/p1
grep -v ^MXE $OUTDIR/$PROJECTTITLE.alltype.exoninfor.txt |cut -f8- >$OUTDIR/p2
grep  ^MXE $OUTDIR/$PROJECTTITLE.alltype.exoninfor.txt  >$OUTDIR/p0
perl -nle '{print "\\t\\t\$_"}' $OUTDIR/p2 >$OUTDIR/pp2

cp $OUTDIR/$PROJECTTITLE.alltype.exoninfor.txt $OUTDIR/$PROJECTTITLE.alltype.exoninfor0.txt

paste $OUTDIR/p1 $OUTDIR/pp2 |cat $OUTDIR/p0 - > $OUTDIR/$PROJECTTITLE.alltype.exoninfor.txt
rm -rf $OUTDIR/p1 $OUTDIR/p2 $OUTDIR/p0 $OUTDIR/pp2
echo Completed!
);


print $cmmd ;
system($cmmd) ;

