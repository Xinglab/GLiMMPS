#!/usr/bin/perl -w

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


open ( IN, "<$IDFILE");
my @ids ;
while ( my $line = <IN>  ){
  chomp($line) ;
  push @ids, $line ;
}
close (IN) ;


mkdir "tempconfig", 0777 unless -d "tempconfig";

$SAMPATHstr= $SAMPATH;
$SAMPATHstr=~s/\//\\\//g;
#print $SAMPATHstr ;


print "Write temporary configuration files to folder tempconfig/...\n";

open (OUT , ">submit_ASreadcounts");

foreach my $BASE1 (@ids) {

  #print $BASE1,"\n" ;
  #print "Write temporary configuration file: tempconfig/config.$BASE1.junctions.txt\n";
  `cp $SCRIPTDIR/config.single.rMATS.txt tempconfig/config.$BASE1.junctions.txt`;
  `perl -pi -e 's/BASE1/$BASE1/g' tempconfig/config.$BASE1.junctions.txt`;
  `perl -pi -e 's/sampath/$SAMPATHstr/g' tempconfig/config.$BASE1.junctions.txt`;
  `perl -pi -e 's/unique\.sam/$SAMFILE/g'  tempconfig/config.$BASE1.junctions.txt`;
  `perl -pi -e 's/PROJECTTITLE/$PROJECTTITLE/g'  tempconfig/config.$BASE1.junctions.txt`;
  `perl -pi -e 's/JUNCTIONLENGTH/$JUNCTIONLENGTH/g'  tempconfig/config.$BASE1.junctions.txt`;
  `perl -pi -e 's/READLENGTH/$READLENGTH/g'  tempconfig/config.$BASE1.junctions.txt`;


  #$cmmd = "perl -p -i -e 's/PROJECTTITLE/$PROJECTTITLE/g'  tempconfig/config.$BASE1.junctions.txt" ;
  #print $cmmd ;
  #system($cmmd) ;
  
  
  print OUT "
python $SCRIPTDIR/rMATS.processsUnique.singlesam.py tempconfig/config.$BASE1.junctions.txt
" ;
###########
## To make a qsub script to run on a cluster, write like this:
## qsub ~/pipelines/sge_Arg1 /opt/Python-2.7.2/bin/python $SCRIPTDIR/rMATS.processsUnique.singlesam.py tempconfig/config.$BASE1.junctions.txt

}
close(OUT) ;

print "Write a tmp script file : submit_ASreadcounts\n\n";


system("chmod +x submit_ASreadcounts");

print "Running submit_ASreadcounts...\n\n";

#system("./submit_ASreadcounts");
