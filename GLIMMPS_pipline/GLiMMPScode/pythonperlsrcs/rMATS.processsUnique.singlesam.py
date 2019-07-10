#
## It takes two samples for the comparison and detection of AS events
## Each sample contains a sam file
#
## MATS input files are generated for each AS events in two different versions
#

### import necessary libraries
import re,os,sys,logging,time,datetime,commands;
import scipy,math;
from scipy import stats;

### checking out the number of arguments
if (len(sys.argv)<2):
  print('Not enough arguments!!');
  print ('It takes one argument');
  print ('Usage:\n\tpython ProgramName.py configFile');
  print ('Example\n\tpython ProgramName.py config.231OE');
  sys.exit();

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;


### setting up the logging format
#logging.basicConfig(level=logging.DEBUG,
#                    format='%(asctime)s %(message)s',
#                    filename='log.processUniqSam.'+  str(datetime.datetime.now()),
#                    filemode='w')

##### Getting Start Time ######
#logging.debug('Start the program with [%s]\n', listToString(sys.argv));
#startTime = time.time();

### setting file pointers
cFile = open(sys.argv[1]); ## config file
#oFile = open('commands.txt', 'a'); ## file that will contain list of commands excuted here

### explain what commands.txt file is for
#oFile.write("#\n### List of commands excuted in the pipeline ###\n#\n");


#### global variables  ################
#### default configuration values #####
readLength=50;
junctionLength=84;
SE='SE';
MXE='MXE';
A5SS='A5SS';
A3SS='A3SS';
AFE='AFE';
ALE='ALE';
RI='RI';
experiment = 'experiment';
base_1 = 'base_1';
dataType = 'paired'; ## either single or paired
samDir = '.';
input_1='.';
outDir = 'myOutput';
email = 'yourID@your.domain';
#######################################

#logging.debug("======================= configuration ============================");
#### assigning proper parameters
#### according to the configuration file
for line in cFile: ## for each line
  if len(line.strip())>0 and line.strip()[0]!='#' : ### it is not the empty line or comments
    ele=line.strip().split('=');
    param = ele[0].replace(' ','')[1:];
    value = ele[1].replace(' ','');
    if param == 'SE': ## Skipped exon file
      SE = value;
#      logging.debug("SE is %s" % SE);
    elif param == 'MXE': ## Mutually exclusive exon
      MXE = value;
#      logging.debug("MXE is %s" % MXE);
    elif param == 'A5SS': ## Alternative 5 prime splice site
      A5SS = value;
#      logging.debug("A5SS is %s" % A5SS);
    elif param == 'A3SS': ## Alternative 3 prime splice site
      A3SS = value;
#      logging.debug("A3SS is %s" % A3SS);
    elif param == 'AFE': ## Alternative first exon
      AFE = value;
#      logging.debug("AFE is %s" % AFE);
    elif param == 'ALE': ## Alternative last exon
      ALE = value;
#      logging.debug("ALE is %s" % ALE);
    elif param == 'RI': ## Retained intron
      RI = value;
#      logging.debug("RI is %s" % RI);
    elif param == 'experiment': ### experiment name. e.g., MB231OE
      experiment = value;
#      logging.debug("experiment or prefix is %s" % experiment);
    elif param == 'base_1': ### base for sample_1. e.g., ESRP
      base_1 = value;
#      logging.debug("base for sample_1 is %s" % base_1);
    elif param == 'base_2': ### base for sample_2. e.g., EV
      base_2 = value;
#      logging.debug("base for sample_2 is %s" % base_2);
    elif param == 'dataType': ## setting data type here, single or paired
      dataType = value;
#      logging.debug("dataType is %s" % dataType);
    elif param == 'samDir': ## sam file directory
      samDir = value;
      if len(samDir.strip())>0: ## we do have samDir
        samDir = samDir.strip()+'/';
#      logging.debug("samDir is %s" % samDir);
    elif param == 'outDir': ## out file directory
      outDir = value;
#      logging.debug("outDir is %s" % outDir);
    elif param == 'input_1': ## sam file names for sample_1
      input_1 = value;
#      logging.debug("input_1 is %s" % input_1);
#    elif param == 'input_2': ## sam file names for sample_2
#      input_2 = value;
#      logging.debug("input_2 is %s" % input_2);
    elif param == 'email': ### alert email
      email = value;
#      logging.debug("email is %s" % email);
    elif param == 'readLength': ### readLength
      readLength = int(value);
#      logging.debug("readLength is %s" % readLength);
    elif param == 'junctionLength': ### junctionLength
      junctionLength = int(value);
#      logging.debug("junctionLength is %s" % junctionLength);

#
commands.getstatusoutput('mkdir '+outDir);
#

## setting up the logging format
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outDir+'/'+'log.processUniqSam.'+  str(datetime.datetime.now()),
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

sample_1 = input_1.split(','); ## could be multiple sam files
#sample_2 = input_2.split(','); ## could be multiple sam files
#
numRep_1 = len(sample_1);
#numRep_2 = len(sample_2);
#
ejLength = junctionLength-readLength+1; ## effective junction length
#
logging.debug("Effective junction length is: %d" % ejLength);
#
logging.debug("There are %d replicates in sample_1: %s" % (numRep_1,base_1));
#logging.debug("There are %d replicates in sample_2: %s" % (numRep_2,base_2));
#
logging.debug("==================================================================");
#
CT1,CT2,CT3=0,1,2; ### count type 1,2, or 3
S1,S2=0,1; ## sample 1 or sample 2
I,S=0,1; ## inclusion isoform or skipping form
#
#
######### more functions #############
#

def logline():
  logging.debug("==================================================================");
## end of logline function

def alertEmail(msg):
  cmd = 'echo "'+msg+'" | mail -s "'+experiment+': processingSam alert" ' + email;
  logging.debug("alert email: %s" % cmd);
  #oFile.write('### sending alert email\n'+cmd+'\n#\n');
  #oFile.flush();
  commands.getstatusoutput(cmd);
### end of alertEmail function

def infoEmail(msg):
  cmd = 'echo "'+msg+'" | mail -s "'+experiment+': processingSam information" ' + email;
  logging.debug("info email: %s" % cmd);
  #oFile.write('### sending information email\n'+cmd+'\n#\n');
  #oFile.flush();
  commands.getstatusoutput(cmd);
### end of infoEmail function

#Input:
#insertsize_inclusion: a numerical variable for the paired-end read inclusion form insert size
#insertsize_skipping: a numerical variable for the paired-end read skipping form insert size
#insertsize_mean: a numerical variable for the mean of the insert size
#insertsize_var: a numerical variable for the variance of the insert size
#Output:
#A 2-element vector for the [fraction of the inclusion form, fraction of the skipping form]
def PE_fraction(insertsize_inclusion, insertsize_skipping, insertsize_mean, insertsize_var):
  p1=stats.norm.cdf(insertsize_inclusion+0.5,insertsize_mean,math.sqrt(insertsize_var));
  p1=p1-stats.norm.cdf(insertsize_inclusion-0.5,insertsize_mean,math.sqrt(insertsize_var));
  p2=stats.norm.cdf(insertsize_skipping+0.5,insertsize_mean,math.sqrt(insertsize_var));
  p2=p2-stats.norm.cdf(insertsize_skipping-0.5,insertsize_mean,math.sqrt(insertsize_var));
  return([p1/(p1+p2),p2/(p1+p2)]);

def getInitialCounts(): ## getting initial counts for each AS event
  rValue = [[[[],[]],[[],[]]],[[[],[]],[[],[]]],[[[],[]],[[],[]]]]; ## count type 1, 2, and 3
  for i in range(0,numRep_1): ## for sample 1
    rValue[CT1][S1][I].append(0);
    rValue[CT1][S1][S].append(0);
    rValue[CT2][S1][I].append(0);
    rValue[CT2][S1][S].append(0);
    rValue[CT3][S1][I].append(0);
    rValue[CT3][S1][S].append(0);
#  for i in range(0,numRep_2): ## for sample 2
#    rValue[CT1][S2][I].append(0);
#    rValue[CT1][S2][S].append(0);
#    rValue[CT2][S2][I].append(0);
#    rValue[CT2][S2][S].append(0);
#    rValue[CT3][S2][I].append(0);
#    rValue[CT3][S2][S].append(0);
  return rValue;
#### end of getInitialCounts()

prefix = experiment;
#
### open AS event files..
#
seFile = open(SE); ## skipped exon event file
mxeFile = open(MXE); ## mxe event file
a5ssFile = open(A5SS); ## a5ss event file
a3ssFile = open(A3SS); ## a3ss event file
#afeFile = open(AFE); ## afe event file
#aleFile = open(ALE); ## ale event file
riFile = open(RI); ## ri event file
#
### open output files here...
#
JC_seFile = open(outDir+'/JC.'+prefix+'.SE.MATS.input.txt', 'w');
JC_mxeFile = open(outDir+'/JC.'+prefix+'.MXE.MATS.input.txt', 'w');
JC_a5ssFile = open(outDir+'/JC.'+prefix+'.A5SS.MATS.input.txt', 'w');
JC_a3ssFile = open(outDir+'/JC.'+prefix+'.A3SS.MATS.input.txt', 'w');
#JC_afeFile = open(outDir+'/JC.'+prefix+'.AFE.MATS.input.txt', 'w');
#JC_aleFile = open(outDir+'/JC.'+prefix+'.ALE.MATS.input.txt', 'w');
JC_riFile = open(outDir+'/JC.'+prefix+'.RI.MATS.input.txt', 'w');
#
JCEC_seFile = open(outDir+'/JCEC.'+prefix+'.SE.MATS.input.txt', 'w');
JCEC_mxeFile = open(outDir+'/JCEC.'+prefix+'.MXE.MATS.input.txt', 'w');
JCEC_a5ssFile = open(outDir+'/JCEC.'+prefix+'.A5SS.MATS.input.txt', 'w');
JCEC_a3ssFile = open(outDir+'/JCEC.'+prefix+'.A3SS.MATS.input.txt', 'w');
#JCEC_afeFile = open(outDir+'/JCEC.'+prefix+'.AFE.MATS.input.txt', 'w');
#JCEC_aleFile = open(outDir+'/JCEC.'+prefix+'.ALE.MATS.input.txt', 'w');
JCEC_riFile = open(outDir+'/JCEC.'+prefix+'.RI.MATS.input.txt', 'w');
#
#JCECPE_seFile = open(outDir+'/JCECPE.'+prefix+'.SE.rMATS.input.txt', 'w');
#JCECPE_mxeFile = open(outDir+'/JCECPE.'+prefix+'.MXE.rMATS.input.txt', 'w');
#JCECPE_a5ssFile = open(outDir+'/JCECPE.'+prefix+'.A5SS.rMATS.input.txt', 'w');
#JCECPE_a3ssFile = open(outDir+'/JCECPE.'+prefix+'.A3SS.rMATS.input.txt', 'w');
#JCECPE_afeFile = open(outDir+'/JCECPE.'+prefix+'.AFE.rMATS.input.txt', 'w');
#JCECPE_aleFile = open(outDir+'/JCECPE.'+prefix+'.ALE.rMATS.input.txt', 'w');
#JCECPE_riFile = open(outDir+'/JCECPE.'+prefix+'.RI.rMATS.input.txt', 'w');
#
logging.debug("populating AS events dictionary");
#
chunk=1000; ## to speed up the sam file processing
#
se={};mxe={};a5ss={};a3ss={};afe={};ale={};ri={}; ## 7 dictionaries
e_se={};e_mxe={};e_a5ss={};e_a3ss={};e_afe={};e_ale={};e_ri={}; ## exons dictionaries
c_se={};c_mxe={};c_a5ss={};c_a3ss={};c_afe={};c_ale={};c_ri={}; ## count dictionaries
s_se={};s_mxe={};s_a5ss={};s_a3ss={};s_afe={};s_ale={};s_ri={}; ## supple dict. effective length of inclusion or skipping form for CT1,CT2,CT3
#
logline();
#
#### SE ######
#
c=0;     ## count
numSE=0; ## number of SE
numSEDup=0; ## duplicate SE id
line=seFile.readline(); ## skipping header
for line in seFile: ## process skipped exon events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
  uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
  dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord

  tLen=min(tE-tS,junctionLength/2); ## target exon length
  uLen=min(uE-uS,junctionLength/2); ## upstream exon length
  dLen=min(dE-dS,junctionLength/2); ## downstream exon length

  e_se[id] = [tS,tE,uS,uE,dS,dE];
  c_se[id] = getInitialCounts();
  I_0= max(0,tLen+uLen-readLength+1)+max(0,tLen+dLen-readLength+1); ## effective inclusion form length for JC
  S_0= max(0,uLen+dLen-readLength+1); ## effective skipping form length for JC
  I_1= max(tE-tS-readLength+1,0)+I_0; ## effective inclusion form length for JC+reads on target
  S_1= S_0; ## effective skipping form length for JC+reads on target
  I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

  #s_se[id] = [[2*ejLength,ejLength],[max(tE-tS-readLength+1,0)+2*ejLength,ejLength],[I_2,S_2]]; ## effective length for CT1,CT2,CT3
  s_se[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

  group = range(uS/chunk, uE/chunk+1) +  range(tS/chunk, tE/chunk+1) +  range(dS/chunk, dE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in se: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in se[chr]: ## this group is already there
        if id in se[chr][i]: ## not likely but this group already has the id
          numSEDup+=1;
          logging.debug("Duplicate SE ID: %d" % id);
        else: ## new SE ID
          se[chr][i][id] = [tS,tE,uS,uE,dS,dE]; ## skipping event with coords
          numSE+=1;
      else: ## new group to this chromosome
        se[chr][i]={};
        se[chr][i][id] = [tS,tE,uS,uE,dS,dE]; ## skipping event with coords
        numSE+=1;
  else: ## first time accesing this chromosome
    se[chr]={};
    for i in group: ## for each possible group
      se[chr][i]={};
      se[chr][i][id] = [tS,tE,uS,uE,dS,dE]; ## skipping event with coords
      numSE+=1;
logging.debug("Processed %d skipped exon events from input AS event file" % c);
logging.debug("Done populating skipped exon dictionary with %d items. %d have duplicate ids" % (numSE, numSEDup));
logging.debug("There are %d se ids in count dictionary and %d se ids in supple dictionary" % (len(c_se),len(s_se)));
#
logline();
#
#### MXE ####
#
c=0;     ## count
numMXE=0; ## number of MXE
numMXEDup=0; ## duplicate MXE id
line=mxeFile.readline(); ## mxe header
for line in mxeFile: ## process mxe events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4]; ## '+' or '-'
  tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
  sS = int(ele[7]); sE = int(ele[8]); ## second exon coord
  if strand=='-': ## negative strand, switch target exon and second exon
    tS = int(ele[7]); tE = int(ele[8]); ## target exon coord
    sS = int(ele[5]); sE = int(ele[6]); ## second exon coord
  uS = int(ele[9]); uE = int(ele[10]); ## upstream exon coord (samller coord)
  dS = int(ele[11]); dE = int(ele[12]); ## downstream exon coord (bigger coord)

  tLen=min(tE-tS,junctionLength/2); ## target exon length
  sLen=min(sE-sS,junctionLength/2); ## second exon length
  uLen=min(uE-uS,junctionLength/2); ## upstream exon length
  dLen=min(dE-dS,junctionLength/2); ## downstream exon length

  e_mxe[id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## target, second, up, down (This is different from the input file)
  c_mxe[id] = getInitialCounts();
  I_0= max(0,tLen+uLen-readLength+1)+max(0,tLen+dLen-readLength+1); ## effective inclusion form length for JC
  S_0= max(0,sLen+uLen-readLength+1)+max(0,sLen+dLen-readLength+1); ## effective skipping form length for JC
  I_1= max(tE-tS-readLength+1,0)+I_0; ## effective inclusion form length for JC+reads on target
  S_1= max(sE-sS-readLength+1,0)+S_0; ## effective skipping form length for JC+reads on target
  I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

  s_mxe[id] = [[2*ejLength,2*ejLength],[max(tE-tS-readLength+1,0)+2*ejLength,max(sE-sS-readLength+1,0)+2*ejLength],[0,0]]; ## effective length for CT1,CT2,CT3
  s_mxe[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

  group = range(uS/chunk,uE/chunk+1)+range(tS/chunk,tE/chunk+1);
  group = group + range(sS/chunk,sE/chunk+1)+range(dS/chunk,dE/chunk+1);## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in mxe: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in mxe[chr]: ## this group is already there
        if id in mxe[chr][i]: ## not likely but this group already has the id
          numMXEDup+=1;
          logging.debug("Duplicate MXE ID: %d" % id);
        else: ## new MXE ID
          mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## mxe event with coords
          numMXE += 1;
      else: ## new group to this chromosome
        mxe[chr][i]={};
        mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## mxe event with coords
        numMXE += 1;
  else: ## first time accesing this chromosome
    mxe[chr]={};
    for i in group: ## for each possible group
      mxe[chr][i]={};
      mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## mxe event with coords
      numMXE+=1;
logging.debug("Processed %d mxe events from input AS event file" % c);
logging.debug("Done populating mxe event dictionary with %d items. %d have duplicate ids" % (numMXE, numMXEDup));
logging.debug("There are %d mxe ids in count dictionary and %d mxe ids in supple dictionary" % (len(c_mxe),len(s_mxe)));
#
logline();
#
#### A5SS ####
#
c=0;     ## count
numA5SS=0; ## number of A5SS
numA5SSDup=0; ## duplicate A5SS id
line=a5ssFile.readline(); ## skipping header
for line in a5ssFile: ## process a5ss events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
  sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord

  lLen=min(lE-lS,junctionLength/2); ## long exon length
  sLen=min(sE-sS,junctionLength/2); ## short exon length
  fLen=min(fE-fS,junctionLength/2); ## flanking exon length
  aLen=min(lE-lS-(sE-sS),junctionLength/2); ## alternative SS region length

  e_a5ss[id] = [lS,lE,sS,sE,fS,fE];
  c_a5ss[id] = getInitialCounts();
  I_0= max(0,lLen+fLen-readLength+1)+max(0,sLen+aLen-readLength+1); ## effective inclusion form length for JC
  S_0= max(0,sLen+fLen-readLength+1); ## effective skipping form length for JC
  I_1= max(lE-lS-(sE-sS)-readLength+1,0)+I_0; ## effective inclusion form length for JC+reads on target
  S_1= S_0; ## effective skipping form length for JC+reads on target
  I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

  #s_a5ss[id] = [[ejLength,ejLength],[max(lE-lS-sE+sS,0)+ejLength,ejLength],[0,0]]; ## effective length for CT1,CT2,CT3
  s_a5ss[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

  group = range(lS/chunk, lE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in a5ss: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in a5ss[chr]: ## this group is already there
        if id in a5ss[chr][i]: ## not likely but this group already has the id
          numA5SSDup+=1;
          logging.debug("Duplicate A5SS ID: %d" % id);
        else: ## new A5SS ID
          a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a5ss event with coords
          numA5SS+=1;
      else: ## new group to this chromosome
        a5ss[chr][i]={};
        a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a5ss event with coords
        numA5SS+=1;
  else: ## first time accesing this chromosome
    a5ss[chr]={};
    for i in group: ## for each possible group
      a5ss[chr][i]={};
      a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a5ss event with coords
      numA5SS+=1;
logging.debug("Processed %d a5ss events from input AS event file" % c);
logging.debug("Done populating a5ss dictionary with %d items. %d have duplicate ids" % (numA5SS, numA5SSDup));
logging.debug("There are %d a5ss ids in count dictionary and %d a5ss ids in supple dictionary" % (len(c_a5ss),len(s_a5ss)));
#
logline();
#
#### A3SS ####
#
c=0;     ## count
numA3SS=0; ## number of A3SS
numA3SSDup=0; ## duplicate A3SS id
line=a3ssFile.readline(); ## skipping header
for line in a3ssFile: ## process a3ss events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
  sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord

  lLen=min(lE-lS,junctionLength/2); ## long exon length
  sLen=min(sE-sS,junctionLength/2); ## short exon length
  fLen=min(fE-fS,junctionLength/2); ## flanking exon length
  aLen=min(lE-lS-(sE-sS),junctionLength/2); ## alternative SS region length

  e_a3ss[id] = [lS,lE,sS,sE,fS,fE];
  c_a3ss[id] = getInitialCounts();
  I_0= max(0,lLen+fLen-readLength+1)+max(0,sLen+aLen-readLength+1); ## effective inclusion form length for JC
  S_0= max(0,sLen+fLen-readLength+1); ## effective skipping form length for JC
  I_1= max(lE-lS-(sE-sS)-readLength+1,0)+I_0; ## effective inclusion form length for JC+reads on target
  S_1= S_0; ## effective skipping form length for JC+reads on target
  I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

  #s_a3ss[id] = [[ejLength,ejLength],[max(lE-lS-sE+sS,0)+ejLength,ejLength],[0,0]]; ## effective length for CT1,CT2,CT3
  s_a3ss[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

  group = range(lS/chunk, lE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in a3ss: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in a3ss[chr]: ## this group is already there
        if id in a3ss[chr][i]: ## not likely but this group already has the id
          numA3SSDup+=1;
          logging.debug("Duplicate A3SS ID: %d" % id);
        else: ## new A3SS ID
          a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a3ss event with coords
          numA3SS+=1;
      else: ## new group to this chromosome
        a3ss[chr][i]={};
        a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a3ss event with coords
        numA3SS+=1;
  else: ## first time accesing this chromosome
    a3ss[chr]={};
    for i in group: ## for each possible group
      a3ss[chr][i]={};
      a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a3ss event with coords
      numA3SS+=1;
logging.debug("Processed %d a3ss events from input AS event file" % c);
logging.debug("Done populating a3ss dictionary with %d items. %d have duplicate ids" % (numA3SS, numA3SSDup));
logging.debug("There are %d a3ss ids in count dictionary and %d a3ss ids in supple dictionary" % (len(c_a3ss),len(s_a3ss)));
#
logline();
#
#### AFE ####
#
#c=0;     ## count
#numAFE=0; ## number of AFE
#numAFEDup=0; ## duplicate AFE id
#line=afeFile.readline(); ## skipping header
#for line in afeFile: ## process afe events file
#  c+=1;
#  ele = line.strip().split('\t');
#  id = int(ele[0]);
#  chr = ele[3];
#  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
#    chr = 'chr'+chr;
#  strand = ele[4];
#  dS = int(ele[5]); dE = int(ele[6]); ## distal exon coord
#  pS = int(ele[7]); pE = int(ele[8]); ## proximal exon coord
#  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord
#
#  dLen=min(dE-dS,junctionLength/2); ## distal exon length
#  pLen=min(pE-pS,junctionLength/2); ## proximal exon length
#  fLen=min(fE-fS,junctionLength/2); ## flanking exon length
#
#  e_afe[id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
#  c_afe[id] = getInitialCounts();
#  I_0= max(0,dLen+fLen-readLength+1); ## effective inclusion form length for JC
#  S_0= max(0,pLen+fLen-readLength+1); ## effective skipping form length for JC
#  I_1= max(dE-dS-readLength+1,0)+I_0; ## effective inclusion form length for JC+reads on target
#  S_1= max(pE-pS-readLength+1,0)+S_0; ## effective skipping form length for JC+reads on target
#  I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet
#
#  #s_afe[id] = [[ejLength,ejLength],[max(dE-dS-readLength+1,0)+ejLength,max(pE-pS-readLength+1,0)+ejLength],[0,0]]; ## effective length for CT1,CT2,CT3
#  s_afe[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3
#
#  group = range(dS/chunk, dE/chunk+1) + range(pS/chunk, pE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
#  group = list(set(group));  ## remove duplicate groups
#
#  if chr in afe: ## already processed this chromosome
#    for i in group: ## for each possible group
#      if i in afe[chr]: ## this group is already there
#        if id in afe[chr][i]: ## not likely but this group already has the id
#          numAFEDup+=1;
#          logging.debug("Duplicate AFE ID: %d" % id);
#        else: ## new AFE ID
#          afe[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
#          numAFE+=1;
#      else: ## new group to this chromosome
#        afe[chr][i]={};
#        afe[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
#        numAFE+=1;
#  else: ## first time accesing this chromosome
#    afe[chr]={};
#    for i in group: ## for each possible group
#      afe[chr][i]={};
#      afe[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
#      numAFE+=1;
#logging.debug("Processed %d afe events from input AS event file" % c);
#logging.debug("Done populating afe dictionary with %d items. %d have duplicate ids" % (numAFE, numAFEDup));
#logging.debug("There are %d afe ids in count dictionary and %d afe ids in supple dictionary" % (len(c_afe),len(s_afe)));
##
#logline();
##
##### ALE ####
##
#c=0;     ## count
#numALE=0; ## number of ALE
#numALEDup=0; ## duplicate ALE id
#line=aleFile.readline(); ## skipping header
#for line in aleFile: ## process ale events file
#  c+=1;
#  ele = line.strip().split('\t');
#  id = int(ele[0]);
#  chr = ele[3];
#  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
#    chr = 'chr'+chr;
#  strand = ele[4];
#  dS = int(ele[5]); dE = int(ele[6]); ## distal exon coord
#  pS = int(ele[7]); pE = int(ele[8]); ## proximal exon coord
#  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord
#
#  dLen=min(dE-dS,junctionLength/2); ## distal exon length
#  pLen=min(pE-pS,junctionLength/2); ## proximal exon length
#  fLen=min(fE-fS,junctionLength/2); ## flanking exon length
#
#  e_ale[id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
#  c_ale[id] = getInitialCounts();
#  I_0= max(0,pLen+fLen-readLength+1); ## effective inclusion form length for JC
#  S_0= max(0,dLen+fLen-readLength+1); ## effective skipping form length for JC
#  I_1= max(pE-pS-readLength+1,0)+I_0; ## effective inclusion form length for JC+reads on target
#  S_1= max(dE-dS-readLength+1,0)+S_0; ## effective skipping form length for JC+reads on target
#  I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet
#
#  #s_ale[id] = [[ejLength,ejLength],[max(pE-pS-readLength+1,0)+ejLength,max(dE-dS-readLength+1,0)+ejLength],[0,0]]; ## effective length for CT1,CT2,CT3
#  s_ale[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3
#
#  group = range(dS/chunk, dE/chunk+1) + range(pS/chunk, pE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
#  group = list(set(group));  ## remove duplicate groups
#
#  if chr in ale: ## already processed this chromosome
#    for i in group: ## for each possible group
#      if i in ale[chr]: ## this group is already there
#        if id in ale[chr][i]: ## not likely but this group already has the id
#          numALEDup+=1;
#          logging.debug("Duplicate ALE ID: %d" % id);
#        else: ## new ALE ID
#          ale[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
#          numALE+=1;
#      else: ## new group to this chromosome
#        ale[chr][i]={};
#        ale[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
#        numALE+=1;
#  else: ## first time accesing this chromosome
#    ale[chr]={};
#    for i in group: ## for each possible group
#      ale[chr][i]={};
#      ale[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
#      numALE+=1;
#logging.debug("Processed %d ale events from input AS event file" % c);
#logging.debug("Done populating ale dictionary with %d items. %d have duplicate ids" % (numALE, numALEDup));
#logging.debug("There are %d ale ids in count dictionary and %d ale ids in supple dictionary" % (len(c_ale),len(s_ale)));
##
#logline();
#
#### RI ####
#
c=0;     ## count
numRI=0; ## number of RI
numRIDup=0; ## duplicate RI id
line=riFile.readline(); ## skipping header
for line in riFile: ## process ri events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  rS = int(ele[5]); rE = int(ele[6]); ## ri exon coord (including up- and down-stream exons)
  uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
  dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord

  rLen=min(rE-rS,junctionLength/2); ## ri exon length
  uLen=min(uE-uS,junctionLength/2); ## upstream exon length
  dLen=min(dE-dS,junctionLength/2); ## downstream exon length
  riLen=min(rE-rS-(uE-uS)-(dE-dS),junctionLength/2); ## retained exon length

  e_ri[id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
  c_ri[id] = getInitialCounts();
  I_0= max(uLen+riLen-readLength+1,0)+max(dLen+riLen-readLength+1,0); ## effective inclusion form length for JC
  S_0= max(uLen+dLen-readLength+1,0); ## effective skipping form length for JC
  I_1= max(rE-rS-(uE-uS)-(dE-dS)-readLength+1,0)+I_0; ## effective inclusion form length for JC+reads on target
  S_1= S_0; ## effective skipping form length for JC+reads on target
  I_2=0; S_2=0; ## effective length considering paired-end info.. not implemented yet

  #s_ri[id] = [[ejLength,ejLength],[dS-uE+readLength-1,ejLength],[0,0]]; ## effective length for CT1,CT2,CT3
  s_ri[id] = [[I_0,S_0],[I_1,S_1],[I_2,S_2]]; ## effective length for CT1,CT2,CT3

  group = range(rS/chunk, rE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in ri: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in ri[chr]: ## this group is already there
        if id in ri[chr][i]: ## not likely but this group already has the id
          numRIDup+=1;
          logging.debug("Duplicate RI ID: %d" % id);
        else: ## new RI ID
          ri[chr][i][id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
          numRI+=1;
      else: ## new group to this chromosome
        ri[chr][i]={};
        ri[chr][i][id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
        numRI+=1;
  else: ## first time accesing this chromosome
    ri[chr]={};
    for i in group: ## for each possible group
      ri[chr][i]={};
      ri[chr][i][id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
      numRI+=1;
logging.debug("Processed %d ri events from input AS event file" % c);
logging.debug("Done populating ri dictionary with %d items. %d have duplicate ids" % (numRI, numRIDup));
logging.debug("There are %d ri ids in count dictionary and %d ri ids in supple dictionary" % (len(c_ri),len(s_ri)));
#
logline();
#
## sys.exit(0);
#
#
def processSample(sample, sInd): ## call it with processSample(sample_1, S1) something like this

  ### process the given sample ###
  for s1 in sample: ## for each sam file
    rep = sample.index(s1);
    if len(s1.strip())<1: ## incorrect split. a user might accidently put comma at the end of input sam file list
      continue; ### just skip this entry, probably the last one though
    sFile = open(samDir+s1.strip()); ## open sam file
    e1 = {}; ## edge count here
    for line in sFile: ## process each line
      if len(line.strip().split('\t'))<5 or line[0]=='#' or line[0]=='@' : ## blank line or comment
        continue;  ## go to next line
      ele = line.strip().split('\t');
      chr = ele[2];
      if chr[0:3]!='chr':
        chr = 'chr'+chr;
      mc = int(ele[3]); ## 1 base, mapping coordinate
      mString = ele[5]; ## mapping string, 50M or aMbNcM format
      group = mc/chunk; ## group does not change, it's okay to check only one group for a junction read
      if 'D' in mString or 'I' in mString or 'S' in mString or 'H' in mString or 'P' in mString or 'X' in mString or '=' in mString: ## skip
        continue; ## go to next line

      ### check to see if the line is either exonic read or junction read
      split_mString = mString.split('M');
      tor = 0; ## type of read, 0 nothing, 1 exonic read, 2 junction read
      if len(split_mString)==2:
        tor = 1; ############ exonic read ######
        rL = int(split_mString[0]); ## read length specified in this mapping string
        mec = mc+rL-1; ## mapping end coord

        ## SE ###
        if chr in se: ## this chromosome has se event(s)
          if group in se[chr]: ## this group has skipped exon event(s)
            for c in se[chr][group]: ## for each skipped exon event in this group
              if (mc>se[chr][group][c][0] and mec<=se[chr][group][c][1]): ## read on the target
                c_se[c][CT2][sInd][I][rep]+=1;
        ### end of SE ###

        ### MXE ####
        if chr in mxe: ## this chromosome has mxe event(s)
          if group in mxe[chr]: ## this group has mxe event(s)
            for c in mxe[chr][group]: ## for each mxe event in this group
              if (mc>mxe[chr][group][c][0] and mec<=mxe[chr][group][c][1]): ## read on the target exon
                c_mxe[c][CT2][sInd][I][rep]+=1;
              elif (mc>mxe[chr][group][c][2] and mec<=mxe[chr][group][c][3]): ## read on the second exon
                c_mxe[c][CT2][sInd][S][rep]+=1;
        ## end of MXE ###

        ## A5SS ##
        if chr in a5ss: ## this chromosome has a5ss event(s)
          if group in a5ss[chr]: ## this group has a5ss event(s)
            for c in a5ss[chr][group]: ## for each a5ss event in this group

              if a5ss[chr][group][c][4]>a5ss[chr][group][c][1]: ## positive strand
           #     if (mc>a5ss[chr][group][c][0] and mc<=(a5ss[chr][group][c][3]-(rL-junctionLength/2)+1) and mec<=a5ss[chr][group][c][1] and mec>=(a5ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
            #      c_a5ss[c][CT1][sInd][I][rep]+=1;
             #     c_a5ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a5ss[chr][group][c][3] and mec<=a5ss[chr][group][c][1]): ## exon read supporting target
                  c_a5ss[c][CT2][sInd][I][rep]+=1;

              else: ## negative strand
             #   if (mc>a5ss[chr][group][c][0] and mc<=(a5ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec<=a5ss[chr][group][c][1] and mec>=(a5ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
              #    c_a5ss[c][CT1][sInd][I][rep]+=1;
               #   c_a5ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a5ss[chr][group][c][0] and mec<=a5ss[chr][group][c][2]): ## exon read supporting target
                  c_a5ss[c][CT2][sInd][I][rep]+=1;

        ## end of A5SS ###

        ## A3SS ##
        if chr in a3ss: ## this chromosome has a3ss event(s)
          if group in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][group]: ## for each a3ss event in this group

              if a3ss[chr][group][c][4]>a3ss[chr][group][c][1]: ## negative strand
                #if (mc>a3ss[chr][group][c][0] and mc<=(a3ss[chr][group][c][3]-(rL-junctionLength/2)+1) and mec<=a3ss[chr][group][c][1] and mec>=(a3ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                 # c_a3ss[c][CT1][sInd][I][rep]+=1;
                 # c_a3ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a3ss[chr][group][c][3] and mec<=a3ss[chr][group][c][1]): ## exon read supporting target
                  c_a3ss[c][CT2][sInd][I][rep]+=1;

              else: ## positive strand
                #if (mc>a3ss[chr][group][c][0] and mc<=(a3ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec<=a3ss[chr][group][c][1] and mec>=(a3ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                 # c_a3ss[c][CT1][sInd][I][rep]+=1;
                 # c_a3ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a3ss[chr][group][c][0] and mec<=a3ss[chr][group][c][2]): ## exon read supporting target
                  c_a3ss[c][CT2][sInd][I][rep]+=1;

        ## end of A3SS ###

#        ## AFE ##
#        if chr in afe: ## this chromosome has afe event(s)
#          if group in afe[chr]: ## this group has afe event(s)
#            for c in afe[chr][group]: ## for each afe event in this group
#              if afe[chr][group][c][4]>afe[chr][group][c][1]: ## positive strand
#                if (mc>afe[chr][group][c][0] and mec<=afe[chr][group][c][1]): ## genome read supporting inclusion (first exon)
#                  c_afe[c][CT2][sInd][I][rep]+=1;
#                elif (mc>afe[chr][group][c][2] and mec<=afe[chr][group][c][3]): ## genome read supporting skipping (second exon)
#                  c_afe[c][CT2][sInd][S][rep]+=1;
#              else: ## negative strand (no difference for genome reads)
#                if (mc>afe[chr][group][c][0] and mec<=afe[chr][group][c][1]): ## genome read supporting inclusion (first exon)
#                  c_afe[c][CT2][sInd][I][rep]+=1;
#                elif (mc>afe[chr][group][c][2] and mec<=afe[chr][group][c][3]): ## genome read supporting skipping (second exon)
#                  c_afe[c][CT2][sInd][S][rep]+=1;
#        ## end of AFE ###
#
#        ## ALE ##
#        if chr in ale: ## this chromosome has ale event(s)
#          if group in ale[chr]: ## this group has ale event(s)
#            for c in ale[chr][group]: ## for each ale event in this group
#              if ale[chr][group][c][4]>ale[chr][group][c][1]: ## negative strand
#                if (mc>ale[chr][group][c][2] and mec<=ale[chr][group][c][3]): ## genome read supporting inclusion (first exon)
#                  c_ale[c][CT2][sInd][I][rep]+=1;
#                elif (mc>ale[chr][group][c][0] and mec<=ale[chr][group][c][1]): ## genome read supporting skipping (second exon)
#                  c_ale[c][CT2][sInd][S][rep]+=1;
#              else: ## positive strand (no difference for genome reads)
#                if (mc>ale[chr][group][c][2] and mec<=ale[chr][group][c][3]): ## genome read supporting inclusion (first exon)
#                  c_ale[c][CT2][sInd][I][rep]+=1;
#                elif (mc>ale[chr][group][c][0] and mec<=ale[chr][group][c][1]): ## genome read supporting skipping (second exon)
#                  c_ale[c][CT2][sInd][S][rep]+=1;
#        ## end of ALE ###

        ## RI ##
        if chr in ri: ## this chromosome has ale event(s)
          if group in ri[chr]: ## this group has ri event(s)
            for c in ri[chr][group]: ## for each ri event in this group, strand does not matter for ri events
              #if (mc>ri[chr][group][c][2] and mc<=ri[chr][group][c][4] and mec>ri[chr][group][c][3] and mec<=ri[chr][group][c][5]): ## genome read supporting retained area
              if (mc>ri[chr][group][c][0] and mc<=(ri[chr][group][c][3]-(rL-junctionLength/2)+1) and mec<=ri[chr][group][c][4] and mec>=(ri[chr][group][c][3]+(rL-junctionLength/2))) or (mc>ri[chr][group][c][3] and mc<=(ri[chr][group][c][4]-(rL-junctionLength/2)+1) and mec<=ri[chr][group][c][5] and mec>=(ri[chr][group][c][4]+(rL-junctionLength/2))): ## multi-exon read supporting target
                c_ri[c][CT1][sInd][I][rep]+=1;
                c_ri[c][CT2][sInd][I][rep]+=1;
              if (mc>ri[chr][group][c][3] and mec<=ri[chr][group][c][4]): ## exon read supporting target
                c_ri[c][CT2][sInd][I][rep]+=1;
        ## end of RI ###


      elif len(split_mString)==3: ###### junction read ###########
        tor = 2; ## junction read
        jS = mc+int(split_mString[0])-1; ## 1-base
        jE = mc+ int(split_mString[0])+ int(split_mString[1].split('N')[0])  -1; ## 0-base
        key = chr+'_'+str(jS)+'_'+str(jE)+'_0';
        if key in e1: ## exist!
          e1[key] = e1[key]+1;
        else: ## new junction
          e1[key] = 1;

        ## SE ###
        if chr in se: ## this chromosome has se event(s)
          if group in se[chr]: ## this group has skipped exon event(s)
            for c in se[chr][group]: ## for each skipped exon event in this group, examine if the given junction is part of it
              if (jS==se[chr][group][c][3] and jE==se[chr][group][c][0]) or (jS==se[chr][group][c][1] and jE==se[chr][group][c][4]): ## IJC
                c_se[c][CT1][sInd][I][rep]+=1;
                c_se[c][CT2][sInd][I][rep]+=1;
              elif jS==se[chr][group][c][3] and jE==se[chr][group][c][4]: ## SJC
                c_se[c][CT1][sInd][S][rep]+=1;
                c_se[c][CT2][sInd][S][rep]+=1;
        ### end of SE ###

        ## MXE ###
        if chr in mxe: ## this chromosome has mxe event(s)
          if group in mxe[chr]: ## this group has mxe event(s)
            for c in mxe[chr][group]: ## for each mxe event in this group, examine if the given junction is part of it
              if (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][0]) or (jS==mxe[chr][group][c][1] and jE==mxe[chr][group][c][6]): ## IJC
                c_mxe[c][CT1][sInd][I][rep]+=1;
                c_mxe[c][CT2][sInd][I][rep]+=1;
              elif (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][2]) or (jS==mxe[chr][group][c][3] and jE==mxe[chr][group][c][6]): ## SJC
                c_mxe[c][CT1][sInd][S][rep]+=1;
                c_mxe[c][CT2][sInd][S][rep]+=1;
        ### end of MXE ###

        ## A5SS ###
        if chr in a5ss: ## this chromosome has a5ss event(s)
          if group in a5ss[chr]: ## this group has a5ss event(s)
            for c in a5ss[chr][group]: ## for each a5ss event in this group, examine if the given junction is part of it
              if a5ss[chr][group][c][4]>a5ss[chr][group][c][1]: ## positive strand
                if jS==a5ss[chr][group][c][1] and jE==a5ss[chr][group][c][4]: ## IJC
                  c_a5ss[c][CT1][sInd][I][rep]+=1;
                  c_a5ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a5ss[chr][group][c][3] and jE==a5ss[chr][group][c][4]: ## SJC
                  c_a5ss[c][CT1][sInd][S][rep]+=1;
                  c_a5ss[c][CT2][sInd][S][rep]+=1;
              else: ## negative strand
                if jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][0]: ## IJC
                  c_a5ss[c][CT1][sInd][I][rep]+=1;
                  c_a5ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][2]: ## SJC
                  c_a5ss[c][CT1][sInd][S][rep]+=1;
                  c_a5ss[c][CT2][sInd][S][rep]+=1;
        ### end of A5SS ###

        ## A3SS ###
        if chr in a3ss: ## this chromosome has a3ss event(s)
          if group in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][group]: ## for each a3ss event in this group, examine if the given junction is part of it
              if a3ss[chr][group][c][4]>a3ss[chr][group][c][1]: ## negative strand
                if jS==a3ss[chr][group][c][1] and jE==a3ss[chr][group][c][4]: ## IJC
                  c_a3ss[c][CT1][sInd][I][rep]+=1;
                  c_a3ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a3ss[chr][group][c][3] and jE==a3ss[chr][group][c][4]: ## SJC
                  c_a3ss[c][CT1][sInd][S][rep]+=1;
                  c_a3ss[c][CT2][sInd][S][rep]+=1;
              else: ## positive strand
                if jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][0]: ## IJC
                  c_a3ss[c][CT1][sInd][I][rep]+=1;
                  c_a3ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][2]: ## SJC
                  c_a3ss[c][CT1][sInd][S][rep]+=1;
                  c_a3ss[c][CT2][sInd][S][rep]+=1;
        ### end of A3SS ###

#        ## AFE ###
#        if chr in afe: ## this chromosome has afe event(s)
#          if group in afe[chr]: ## this group has afe event(s)
#            for c in afe[chr][group]: ## for each afe event in this group, examine if the given junction is part of it
#              if afe[chr][group][c][4]>afe[chr][group][c][1]: ## positive strand
#                if jS==afe[chr][group][c][3] and jE==afe[chr][group][c][4]: ## IJC
#                  c_afe[c][CT1][sInd][I][rep]+=1;
#                  c_afe[c][CT2][sInd][I][rep]+=1;
#                elif jS==afe[chr][group][c][1] and jE==afe[chr][group][c][4]: ## SJC
#                  c_afe[c][CT1][sInd][S][rep]+=1;
#                  c_afe[c][CT2][sInd][S][rep]+=1;
#              else: ## negative strand
#                if jS==afe[chr][group][c][5] and jE==afe[chr][group][c][2]: ## IJC
#                  c_afe[c][CT1][sInd][I][rep]+=1;
#                  c_afe[c][CT2][sInd][I][rep]+=1;
#                elif jS==afe[chr][group][c][5] and jE==afe[chr][group][c][0]: ## SJC
#                  c_afe[c][CT1][sInd][S][rep]+=1;
#                  c_afe[c][CT2][sInd][S][rep]+=1;
#        ### end of AFE ###
#
#        ## ALE ###
#        if chr in ale: ## this chromosome has ale event(s)
#          if group in ale[chr]: ## this group has ale event(s)
#            for c in ale[chr][group]: ## for each ale event in this group, examine if the given junction is part of it
#              if ale[chr][group][c][4]>ale[chr][group][c][1]: ## negative strand
#                if jS==ale[chr][group][c][3] and jE==ale[chr][group][c][4]: ## IJC
#                  c_ale[c][CT1][sInd][I][rep]+=1;
#                  c_ale[c][CT2][sInd][I][rep]+=1;
#                elif jS==ale[chr][group][c][1] and jE==ale[chr][group][c][4]: ## SJC
#                  c_ale[c][CT1][sInd][S][rep]+=1;
#                  c_ale[c][CT2][sInd][S][rep]+=1;
#              else: ## positive strand
#                if jS==ale[chr][group][c][5] and jE==ale[chr][group][c][2]: ## IJC
#                  c_ale[c][CT1][sInd][I][rep]+=1;
#                  c_ale[c][CT2][sInd][I][rep]+=1;
#                elif jS==ale[chr][group][c][5] and jE==ale[chr][group][c][0]: ## SJC
#                  c_ale[c][CT1][sInd][S][rep]+=1;
#                  c_ale[c][CT2][sInd][S][rep]+=1;
#        ### end of ALE ###

        ## RI ###
        if chr in ri: ## this chromosome has ri event(s)
          if group in ri[chr]: ## this group has ri event(s)
            for c in ri[chr][group]: ## for each ri event in this group, examine if the given junction is part of it
              if jS==ri[chr][group][c][3] and jE==ri[chr][group][c][4]: ## SJC
                c_ri[c][CT1][sInd][S][rep]+=1;
                c_ri[c][CT2][sInd][S][rep]+=1;
        ### end of RI ###


      else: ## it is not exonic nor junction read. proceed to the next line
        continue;

    logging.debug("Done populating edeg count for %s with %d junctions" % (s1.strip(), len(e1)));
    sFile.close();
    edgeFile = open(samDir+s1.strip()+'.edgeCount', 'w');
    for k in e1:
      edgeFile.write(k+'\t'+str(e1[k])+'\n');
    edgeFile.close();
##### end of processSample ######

if dataType=="skipThisForNow": ## call paired..
  logging.debug("Start processing sample_1: %s" % base_1);
  processSample_PE(sample_1, S1, 72, 60);
  logging.debug("Done processing %s" % base_1);
#  logging.debug("Start processing sample_2: %s" % base_2);
#  processSample_PE(sample_2, S2, 75, 62);
#  logging.debug("Done processing %s" % base_2);
else:
  logging.debug("Start processing sample_1: %s" % base_1);
  processSample(sample_1, S1);
  logging.debug("Done processing %s" % base_1);
#  logging.debug("Start processing sample_2: %s" % base_2);
#  processSample(sample_2, S2);
#  logging.debug("Done processing %s" % base_2);

###
#def writeInputFile(h1,h2,h3,f1,f2,f3,cnt,sup): ## header 1,2,3, file 1,2,3, count dict, supple dict)
def writeInputFile(h1,h2,h3,f1,f2,cnt,sup): ## header 1,2,3, file 1,2,3, count dict, supple dict)
  ## print header first
  f1.write(h1+'\n');
  f2.write(h2+'\n');
#  f3.write(h3+'\n');

  for k in sorted(sup.keys()):
    #f1.write(str(k)+'\t'+str(cnt[k][CT1][S1][I][0])+'\t'+str(cnt[k][CT1][S1][S][0])+'\t'+str(cnt[k][CT1][S2][I][0])+'\t'+str(cnt[k][CT1][S2][S][0])+'\t'+str(sup[k][CT1][0])+'\t'+str(sup[k][CT1][1])+'\n'); ## need to correct this to handle replicates...
    #f2.write(str(k)+'\t'+str(cnt[k][CT2][S1][I][0])+'\t'+str(cnt[k][CT2][S1][S][0])+'\t'+str(cnt[k][CT2][S2][I][0])+'\t'+str(cnt[k][CT2][S2][S][0])+'\t'+str(sup[k][CT2][0])+'\t'+str(sup[k][CT2][1])+'\n'); ## need to correct this to handle replicates...
    f1.write(str(k)+'\t'+','.join(map(str,cnt[k][CT1][S1][I]))+'\t'+','.join(map(str,cnt[k][CT1][S1][S]))+'\t'+str(sup[k][CT1][0])+'\t'+str(sup[k][CT1][1])+'\n'); ## need to correct this to handle replicates...
    f2.write(str(k)+'\t'+','.join(map(str,cnt[k][CT2][S1][I]))+'\t'+','.join(map(str,cnt[k][CT2][S1][S]))+'\t'+str(sup[k][CT2][0])+'\t'+str(sup[k][CT2][1])+'\n'); ## need to correct this to handle replicates...
##### end of writeInputFile function #######

#
logging.debug("Writing out the rMATS input file..");
#
CT1_header = 'ID\tIJC_'+base_1+'\tSJC_'+base_1+'\tIncFormLen\tSkipFormLen';
CT2_header = 'ID\tIC_'+base_1+'\tSC_'+base_1+'\tIncFormLen\tSkipFormLen';
CT3_header = 'ID\tIC_'+base_1+'\tSC_'+base_1+'\tIncFormLen\tSkipFormLen';
#
## SE ##
#writeInputFile(CT1_header,CT2_header,CT3_header,JC_seFile,JCEC_seFile,JCECPE_seFile,c_se,s_se);
writeInputFile(CT1_header,CT2_header,CT3_header,JC_seFile,JCEC_seFile,c_se,s_se);
## MXE ##
#writeInputFile(CT1_header,CT2_header,CT3_header,JC_mxeFile,JCEC_mxeFile,JCECPE_mxeFile,c_mxe,s_mxe);
writeInputFile(CT1_header,CT2_header,CT3_header,JC_mxeFile,JCEC_mxeFile,c_mxe,s_mxe);
## A5SS ##...
#writeInputFile(CT1_header,CT2_header,CT3_header,JC_a5ssFile,JCEC_a5ssFile,JCECPE_a5ssFile,c_a5ss,s_a5ss);
writeInputFile(CT1_header,CT2_header,CT3_header,JC_a5ssFile,JCEC_a5ssFile,c_a5ss,s_a5ss);
## A3SS ##...
#writeInputFile(CT1_header,CT2_header,CT3_header,JC_a3ssFile,JCEC_a3ssFile,JCECPE_a3ssFile,c_a3ss,s_a3ss);
writeInputFile(CT1_header,CT2_header,CT3_header,JC_a3ssFile,JCEC_a3ssFile,c_a3ss,s_a3ss);
### AFE ##
#writeInputFile(CT1_header,CT2_header,CT3_header,JC_afeFile,JCEC_afeFile,JCECPE_afeFile,c_afe,s_afe);
### ALE ##...
#writeInputFile(CT1_header,CT2_header,CT3_header,JC_aleFile,JCEC_aleFile,JCECPE_aleFile,c_ale,s_ale);
## RI ##...
#writeInputFile(CT1_header,CT2_header,CT3_header,JC_riFile,JCEC_riFile,JCECPE_riFile,c_ri,s_ri);
writeInputFile(CT1_header,CT2_header,CT3_header,JC_riFile,JCEC_riFile,c_ri,s_ri);
#
#### close all files here ##########
seFile.close()
mxeFile.close()
a5ssFile.close()
a3ssFile.close()
#afeFile.close()
#aleFile.close()
riFile.close()
#
JC_seFile.close()
JC_mxeFile.close()
JC_a5ssFile.close()
JC_a3ssFile.close()
#JC_afeFile.close()
#JC_aleFile.close()
JC_riFile.close()
#
JCEC_seFile.close()
JCEC_mxeFile.close()
JCEC_a5ssFile.close()
JCEC_a3ssFile.close()
#JCEC_afeFile.close()
#JCEC_aleFile.close()
JCEC_riFile.close()
#
#JCECPE_seFile.close()
#JCECPE_mxeFile.close()
#JCECPE_a5ssFile.close()
#JCECPE_a3ssFile.close()
#JCECPE_afeFile.close()
#JCECPE_aleFile.close()
#JCECPE_riFile.close()
#
######################################
#
infoEmail("Processing unique sam file is finished with no error!");
#
#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
#
sys.exit(0);
#
