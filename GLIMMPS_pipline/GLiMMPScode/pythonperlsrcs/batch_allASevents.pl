#!/usr/bin/perl -w
##　For all individuals in the population, obtain the all the possible Alternative Splicing (AS) events in the population based on either the gene annotation or the new splicing site identified from tophat.


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

open ( IN, "<$IDFILE");  ## The file list all Indivividual IDs in the population.
my @ids = <IN> ;

#make directory unless it already exists
mkdir "ASEvents", 0777 unless -d "ASEvents";
mkdir "temp", 0777 unless -d "temp";



#################
# write to a tmp script file 

print "Write a tmp script file : submit_ASevents\n\n";

open (OUT, ">submit_ASevents");

my $str = "python $SCRIPTDIR/processGTF.SAMs.py $ANNOTATIONFILE $PROJECTTITLE ";
print OUT $str ; 
my $count =0 ;
#print "All individuals:\n";
foreach my $id (@ids) {
chomp($id) ;
#print $id,"\n" ;

$count++ ;
if ($count>1) {print OUT "," ; } 
print OUT "$SAMPATH/$id/$SAMFILE" ;

}
print OUT " temp" ;
close(OUT) ;

system("chmod +x submit_ASevents");
#system("./submit_ASevents");

