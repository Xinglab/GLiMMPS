#
## this program processes  GTF file and SAM files to get both known and novel AS events
#

### import necessary libraries
import re,os,sys,logging,time,datetime;
#
myVer= "2.5";
#
### checking out the number of arguments
if (len(sys.argv)<4): 
  print('Not enough arguments!!');
  print ('It takes at least 3 arguments.');
  print ('Usage:\n\tpython ProgramName.py gtfFile outputPrefix SAMfiles [logFolder]');
  print ('Example\n\tpython ProgramName.py AceView.ncbi_37.gtf fromAceView sample1.sam,sample2.sam,... temp');
  sys.exit();

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;

def uniq(inlist):
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

### setting up the logging format 
logFolder = '';
if len(sys.argv)>4: ## we have a log folder
  logFolder = sys.argv[4].strip()+'/';
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=logFolder+'log.process.GTF.SAMs.'+myVer+'.'+ str(datetime.datetime.now()),
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();
OUTDIR ="ASEvents";
###
iFile = open(sys.argv[1]); ## input gtf file
oFile_3 = open(OUTDIR+'/'+sys.argv[2]+'.A3SS.txt', 'w'); ## alt-3 SS output file
oFile_5 = open(OUTDIR+'/'+sys.argv[2]+'.A5SS.txt', 'w'); ## alt-5 SS output file
oFile_ce = open(OUTDIR+'/'+sys.argv[2]+'.SE.txt', 'w'); ## skipped exon (cassette exon) output file
oFile_mxe = open(OUTDIR+'/'+sys.argv[2]+'.MXE.txt', 'w'); ## Mutually exclusive exon output file
oFile_afe = open(OUTDIR+'/'+sys.argv[2]+'.AFE.txt', 'w'); ## alt first exon output file
oFile_ale = open(OUTDIR+'/'+sys.argv[2]+'.ALE.txt', 'w'); ## alt last exon output file
oFile_ri = open(OUTDIR+'/'+sys.argv[2]+'.RI.txt', 'w'); ## retained intron output file
#
samFiles = sys.argv[3].split(','); ## coulbe be multiple sam files
#


### write out header
ceHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE";
oFile_ce.write(ceHeader+'\n');

mxeHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE";
oFile_mxe.write(mxeHeader+'\n');
#oFile_mxe_filtered.write(mxeHeader+'\n');

altSSHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE";
oFile_3.write(altSSHeader+'\n');
oFile_5.write(altSSHeader+'\n');

altFLHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\tdistalExonStart_0base\tdistalExonEnd\tproximalES\tproximalEE\tflankingES\tflankingEE"
oFile_afe.write(altFLHeader+'\n');
oFile_ale.write(altFLHeader+'\n');

riHeader = "ID\tGeneID\tgeneSymbol\tchr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE";
oFile_ri.write(riHeader+'\n');

c=0;

chunk=1000;
geneGroup={}; ## genes in a group

genes = {}; 
supple = {}; 
cds={}; ## coding region

for line in iFile: ## for each line
  if ( line.startswith('##')): continue; ## skip the comment lines
  ele = line.strip().split('\t');
  chr = ele[0];
  type = ele[2]; ## exon, intron, CDS, start_codon, stop_codon..
  sC = ele[3]; ## start coord, 1-base
  eC = ele[4]; ## end coord, 1-base
  group = range(int(sC)/chunk, int(eC)/chunk + 1); ## groups this line could belong to
  group = list(set(group));  ## remove duplicate groups
  strand = ele[6]; 
  desc = ele[8].split(';')
  gID = desc[0].strip().split(' ');
  txID = desc[1].strip().split(' ');

  if gID[0]!='gene_id' or txID[0] != 'transcript_id': ## wrong one..
    logging.debug("gID or txID is not correct: %s, %s" % (gID, txID));
    logging.debug("Incorrect description: %s" % ele);
    continue; ## process next line

  for i in group: ## for each possible group
    if i in geneGroup: ## this group already exist
      geneGroup[i].append(gID[1]); ## duplicate geneIDs will get removed after the outer for loop
    else: ## first time accesing this group
      geneGroup[i] = [gID[1]];

  if type=='exon':  ## process exon
    if gID[1] in genes: # already processed this gID
      if txID[1] in genes[gID[1]]: ## this transcript is added already
        genes[gID[1]][txID[1]].append([int(sC), int(eC)]); ## add exon to the existing Tx
      else: ## first time processing this Tx
        genes[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
    else:  ## new gene ID
      genes[gID[1]] = {};
      genes[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
      supple[gID[1]] = [gID[1], chr, strand]; ## geneID, chromosom and strand
  if type=='CDS': ## coding region
    if gID[1] in cds: # already processed this gID
      if txID[1] in cds[gID[1]]: ## this transcript is added already
        cds[gID[1]][txID[1]].append([int(sC), int(eC)]); ## add CDS to the existing Tx
      else: ## first time processing this Tx
        cds[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first CDS
    else:  ## new gene ID
      cds[gID[1]] = {};
      cds[gID[1]][txID[1]] = [[int(sC), int(eC)]]; ## add first exon
logging.debug("Done populating genes and cds dictionaries");
#
for gg in geneGroup: ## for all groups in geneGroup
  geneGroup[gg] = list(set(geneGroup[gg]));
#
## stats
#
logging.debug("======== stats from genes =========");
#
nGene=len(genes); ## number of genes in genes dict
nTx=0; ## number of transcripts
oneTx=0; ## number of one-tx genes
nExon = 0; ## number of exons
oneExon=0; ## number of one-exon transcripts
#
oneTxOneExon=0; ## number of one-tx genes with only one exon
#
for id in genes: ## for each gene
  nTx += len(genes[id]); 
  if len(genes[id])==1:
    oneTx += 1; ## one-transcript gene
  for tx in genes[id]: ## for each tx
    nExon += len(genes[id][tx]);
    if len(genes[id][tx])==1: ## one exon tx
      oneExon += 1;
      if len(genes[id])==1: ## one tx gene
        oneTxOneExon+=1;

logging.debug("There are %d distinct gene ID in the gtf file" % nGene);
logging.debug("There are %d distinct transcript ID in the gtf file" % nTx);
logging.debug("There are %d one-transcript genes in the gtf file" % oneTx);
logging.debug("There are %d exons in the gtf file" % nExon);
logging.debug("There are %d one-exon transcripts in the gtf file" % oneExon);
logging.debug("There are %d one-transcript genes with only one exon in the transcript" % oneTxOneExon);
logging.debug("Average number of transcripts per gene is %f" % (float(nTx)/nGene));
logging.debug("Average number of exons per transcript is %f" % (float(nExon)/nTx));
logging.debug("Average number of exons per transcript excluding one-exon tx is %f" % (float(nExon-oneExon)/(nTx-oneExon)));
#
logging.debug("======== stats from cds =========");
#
nGene=len(cds); ## number of genes in cds dict
nTx=0; ## number of transcripts
oneTx=0; ## number of one-tx genes
nExon = 0; ## number of exons
oneExon=0; ## number of one-exon transcripts
#
oneTxOneExon=0; ## number of one-tx genes with only one exon
#
for id in cds: ## for each gene
  nTx += len(cds[id]);
  if len(cds[id])==1:
    oneTx += 1; ## one-transcript gene
  for tx in cds[id]: ## for each tx
    nExon += len(cds[id][tx]);
    if len(cds[id][tx])==1: ## one exon tx
      oneExon += 1;
      if len(cds[id])==1: ## one tx gene
        oneTxOneExon+=1;

logging.debug("There are %d distinct gene ID in the gtf file" % nGene);
logging.debug("There are %d distinct transcript ID in the gtf file" % nTx);
logging.debug("There are %d one-transcript genes in the gtf file" % oneTx);
logging.debug("There are %d exons in the gtf file" % nExon);
logging.debug("There are %d one-exon transcripts in the gtf file" % oneExon);
logging.debug("There are %d one-transcript genes with only one cds in the transcript" % oneTxOneExon);
logging.debug("Average number of transcripts per gene is %f" % (float(nTx)/nGene));
logging.debug("Average number of cds per transcript is %f" % (float(nExon)/nTx));
logging.debug("Average number of cds per transcript excluding one-cds tx is %f" % (float(nExon-oneExon)/(nTx-oneExon)));
#
#
logging.debug("======== stats from geneGroup =========");
#
tgi = 0;## total geneIDs
## 
for gg in geneGroup: 
  tgi += len(geneGroup[gg]);
logging.debug("There are total of %d groups and %d genes in geneGroup" % (len(geneGroup), tgi));
logging.debug("The average number of genes in each group is %f" % (float(tgi)/len(geneGroup))); 
#
logging.debug("==========================================\n");
#
#
### sort transcripts ###
#
for gID in genes:
  for tID in genes[gID]: ## sort each transcript
    if len(genes[gID][tID])==1: ## only one exon, skip it
      continue; ## next transcript..
    genes[gID][tID] = sorted(genes[gID][tID]); ## sort each transcript
#
#
### now process SAM files to add novel transcripts constructed from novel junctions 
### novel junction connects the existing exons (within a gene) but the junction is not defined in the gtf
#
samIndex=0;
for s1 in samFiles: ## for each samFiles
  if len(s1.strip())<1: ## incorrect split. a user might accidently put comma at the end of input sam file list
    continue; ### just skip this entry, probably the last one though
  samIndex+=1; ## for the novel tx index
  logging.debug("processing %s" % s1.strip());
  sFile = open(s1.strip()); ## open sam file
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
      #tor = 1; ## exonic read ##
      continue; ## go to next line
    elif len(split_mString)==3: ###### junction read ###########
      tor = 2; ## junction read
      jS = mc+int(split_mString[0])-1; ## 1-base
      jE = mc+ int(split_mString[0])+ int(split_mString[1].split('N')[0])  -1; ## 0-base
      #key = chr+'_'+str(jS)+'_'+str(jE)+'_0';
  
      ## get a group of genes then examine each gene in the group
      if group not in geneGroup: ## this junction is not from annotated genes
        continue; ## go to the next line
      for cg in geneGroup[group]: ## for each candidate gene
        ntx={}; ## novel transcript for the given candidate gene
        for ctx in genes[cg]: ## for each transcript in the candidate gene
          cexons = genes[cg][ctx]; ## candidate exons
          uInd=-1; dInd=-1; ## upstream exon index and downstream exon index

          for ci in range(0, len(cexons)): ## examine each candidate exon
            if cexons[ci][1]==jS: ## this exon has the same end
              uInd=ci;
            if cexons[ci][0]==(jE+1): ## this exon has the same start
              dInd=ci;

          if uInd>-1 and dInd>-1: ## we have both exons here..
            if dInd-uInd == 1: ## known junction.. examine next gene
              #if jS==18542228 and (jE+1)==18542603:
              #  print cexons;
              #  print ntx, uInd, dInd, samIndex; 
              #  sys.exit();
              ntx={}; ## empty novel transcript dictionary then go to the next cg
              break; ## break for ctx loop.. go and process next cg
            elif dInd-uInd>1: ## novel junction.. make novel transcripts then add to ntx..
              key=str(cexons[uInd][0])+':'+str(cexons[uInd][1])+':'+str(cexons[dInd][0])+':'+str(cexons[dInd][1]); ## novel junction
              nex = [cexons[uInd],cexons[dInd]]; ## novel exons 
              if uInd>0: ## not the first exon
                key = str(cexons[uInd-1][0])+':'+str(cexons[uInd-1][1])+':'+key; ## prev exon + novel junction
                nex = [cexons[uInd-1]]+nex;
              if dInd<len(cexons)-1: ## not the last exon
                key = key + ':'+str(cexons[dInd+1][0])+':'+str(cexons[dInd+1][1]); ## novel junction + next exon
                nex = nex+[cexons[dInd+1]];
              ntx[key] = nex;  ## it's okay to overwrite
              #if jS==18542228 and (jE+1)==18542603:
              #  print ntx, uInd, dInd, samIndex; 
              #print ("Junction: %d - %d" % (jS,jE+1));
              #print ("%s: %s" % (key, ntx[key]));
            else: ### should not be here..
              logging.debug("Check this out..gene: %s, tx: %s, uExon: %s, dExon: %s" % (cg,ctx,cexons[uInd], cexons[dInd]));
        ## end of for ctx in genes[cg]
        for novelT in ntx: ## for all novel transcript
          txName = cg+'.novel_'+str(samIndex)+'_'+novelT;
          #ntIndex += 1;
          genes[cg][txName] = ntx[novelT];
        #if jS==18542228 and (jE+1)==18542603:
        #  print ntx;
        #  print genes;
        #  sys.exit(); 
       
      ## end of for cg in geneGroup[group]

    ## end of elif junction read

  ## end of for line in sFile
  logging.debug("Done processing %s" % s1.strip());
#
#print "=== genes ===";
#for ggg in genes:
#  for ttt in genes[ggg]:
#    print "%s,%s: %s" % (ggg,ttt,genes[ggg][ttt]);
#
def fullyContainedInInternalExon(myExon, myGeneID): ## return true if any transcript has an internal exon fully containing the exon
  for myTxID in genes[myGeneID]:
    for intExon in genes[myGeneID][myTxID][1:-1]: ## for each internal exon
      if intExon[0]<=myExon[0] and intExon[1]>=myExon[1]: ## fully contained
        return True;
  return False;
#
logging.debug("Process each gene from dictionary");
#
numSkippingEvents=0;
numMXEvents=0;
numSSEvents=0;
num3=0;
num5=0;
numAFE=0;
numALE=0;
numRI=0; 

sEvents={};
mxEvents={};
ss3Events={};
ss5Events={};
mxEvents_filtered={}; 
afeEvents={}; ## alternative first exon
aleEvents={}; ## alternative last exon
riEvents={};  ## retained intron 

dupSE=0;
dupMXE=0;
dupSS3=0;
dupSS5=0;
dupRI=0;

filteredMXE=0;
#
for gID in genes:  ## process each gene
  supInfo = supple[gID]; ## supplementary info
  if len(genes[gID])==1: ## only one transcript, alt SS event is imposible
    continue; ## go to the next geneID
  else: ## look for alt SS event
    de={}; ## distinct exons
    for tID in genes[gID]: ## sort each transcript, merge and get distinct exons
      if len(genes[gID][tID])==1: ## only one exon, skip it 
        continue; ## next transcript..
      genes[gID][tID] = sorted(genes[gID][tID]); ## sort each transcript
      for exon in genes[gID][tID]:
        de[exon[0], exon[1]] = 1; ## it's okay to overwrite

    dc={}; ## distinct cds
    if gID in cds: ## this gene has cds 
      for tID in cds[gID]: ## sort each transcript, merge and get distinct cds
        if len(cds[gID][tID])==1: ## only one cds, skip it
          continue; ## next transcript..
        cds[gID][tID] = sorted(cds[gID][tID]); ## sort each transcript
        for exon in cds[gID][tID]: ## for each cds
          dc[exon[0], exon[1]] = 1; ## it's okay to overwrite

    ## now we have sorted transcripts and distinct exons
    ## examine each exon in de to see if the exon is involved in any types of AS events
  
    for ce in de: ## for each exon in distinct exon dictionary
      uf=[]; ## upstream flanking exons
      df=[]; ## downstream flanking exons
      for tID in genes[gID]: ## examine each transcript to see if it contains the given exon
        if [ce[0],ce[1]] in genes[gID][tID]: ## this exon is in the transcript     
          eInd = genes[gID][tID].index([ce[0],ce[1]]);
          if (0<eInd): ## it definitely has upstream flanking exon
            uf.append(genes[gID][tID][eInd-1]);
          if (eInd<len(genes[gID][tID])-1): ## it definitely has downstream flanking exon
            df.append(genes[gID][tID][eInd+1]);

      ## getting uniq upstream flanking exons and downstream flanking exons
      uf = uniq(uf);
      df = uniq(df);

      #### exon skipping (cassette exon) events ###
      for i in range(0,len(uf)): ### going through the upstream flanking exons
        f1 = uf[i]; ## first flanking exon
        for j in range(0,len(df)): ## going through the downstream falnking exons to see if there is skipping junction
          f2 = df[j];
          for tID in genes[gID]: ## examine each transcript, to see if it has an edge connecting f1 and f2
            if len(genes[gID][tID])<2: ## less than two exons, skip it
              continue;

########################## not requiring exact falnking exons ###########
            for i in range(0,len(genes[gID][tID])-1): ## for each exon in the tx
              e_1 = genes[gID][tID][i];
              e_2 = genes[gID][tID][i+1]; 
              if e_1[1]==f1[1] and e_2[0]==f2[0]: ## this tx has an edge connecting f1 and f2 but does not have ce
                key=supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(f1[1])+':'+str(f2[0]-1); ## target exon and skipping junction
                if key in sEvents: ## already have this skipping events
                  dupSE +=1;
                  continue; ## next transcript
                else: ## new key, write it
                  sEvents[key]=1;
                  numSkippingEvents += 1;
                  oFile_ce.write(str(numSkippingEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
###################  

####   requiring exact flanking exons ##################
#            if f1 in genes[gID][tID] and f2 in genes[gID][tID]: ## this transcript has f1 and f2
#              if genes[gID][tID].index(f1)+1 == genes[gID][tID].index(f2): ### there is an edge connecting f1 and f2
#                #key=supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(f2[0]-1)+':'+str(f2[1]);
#                key=supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(f1[1])+':'+str(f2[0]-1); ## target exon and skipping junction
#                if key in sEvents: ## already have this skipping events
#                  dupSE +=1;
#                  continue; ## next transcript
#                else: ## new key, write it
#                  sEvents[key]=1;
#                  numSkippingEvents += 1; 
#                  oFile_ce.write(str(numSkippingEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
##########################################################################

      ##### mutually exclusive events #####
      for i in range(0,len(uf)): ### going through the upstream flanking exons
        f1 = uf[i]; ## first flanking exon
        for j in range(0,len(df)): ## going through the downstream falnking exons to see if there is skipping junction
          f2 = df[j];
          for tID in genes[gID]: ## examine each transcript
            if [ce[0], ce[1]] in genes[gID][tID] or f1 not in genes[gID][tID] or f2 not in genes[gID][tID]: ## ce in, f1 or f2 is not in..do not examine
              continue; ### go to next transcript in genes[gID]
            else: ## this transcript does not have ce and has both f1 and f2, let's take a look
              fromF1 = genes[gID][tID].index(f1);
              fromF2 = genes[gID][tID].index(f2);
              mxe = genes[gID][tID][fromF1+1]; ## candidate mxe
              if (fromF1+1 == fromF2-1) and (mxe[0]>ce[1]): ### this exon is the right one
                goodMXE=True;
#                for txs in genes[gID]: ### search through transcripts again to see if there is ce-mxe edge
#                  if [ce[0],ce[1]] in genes[gID][txs] and mxe in genes[gID][txs]: ## ce and mxe are not MXE
#                    goodMXE = False;
#                    filteredMXE+=1;
#                    fkey = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(mxe[0]-1)+':'+str(mxe[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(f2[0]-1)+':'+str(f2[1]);
#                    if fkey in mxEvents_filtered: ## duplicate, do not write out
#                      pass;
#                    else:
#                      mxEvents_filtered[fkey] = 1;
#                      oFile_mxe_filtered.write(str(filteredMXE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(mxe[0]-1)+'\t'+str(mxe[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');
#                    break; ### not MXE
                if goodMXE: ### it's okay to write out
                  #key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(mxe[0]-1)+':'+str(mxe[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(f2[0]-1)+':'+str(f2[1]);
                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(mxe[0]-1)+':'+str(mxe[1])+':'+str(f1[1])+':'+str(f2[0]-1);
                  if key in mxEvents: ## duplicate, 
                    dupMXE += 1;
                  else:
                    numMXEvents += 1;
                    mxEvents[key] = 1; 
                    oFile_mxe.write(str(numMXEvents)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(mxe[0]-1)+'\t'+str(mxe[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(f2[0]-1)+'\t'+str(f2[1])+'\n');


      #### alt-3 and alt-5 events ###
      for i in range(0,len(uf)-1): ### going through the upstream flanking exons
        e=uf[i]; 
        for j in range(i+1, len(uf)):
          u=uf[j];
          if e[0]==u[0]: ## it is alt SS event, because uf is derived from SET
            #key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(e[0]-1)+':'+str(max(e[1],u[1]))+':'+str(e[0]-1)+':'+str(min(e[1],u[1]));
            key = supInfo[1]+':'+str(min(e[1],u[1]))+':'+str(max(e[1],u[1]))+':'+str(ce[0]-1);
            if supInfo[2]=='+': ## positive strand. alt-5 event
              if key in ss5Events: ## duplicate
                dupSS5 += 1;
              else:
                ss5Events[key]=1;
                num5 += 1;
                oFile_5.write(str(num5)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(e[0]-1)+'\t'+str(max(e[1],u[1]))+'\t'+str(e[0]-1)+'\t'+str(min(e[1],u[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
            else: ## neg strand. alt-3 event
              if key in ss3Events: ## duplicate
                dupSS3 += 1;
              else:
                ss3Events[key]=1;
                num3 += 1;
                oFile_3.write(str(num3)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(e[0]-1)+'\t'+str(max(e[1],u[1]))+'\t'+str(e[0]-1)+'\t'+str(min(e[1],u[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');

      for i in range(0,len(df)-1): ### going through the downstream flanking exons
        e=df[i];
        for j in range(i+1, len(df)):
          d=df[j];
          if e[1]==d[1]: ## it is alt SS event, because uf is derived from SET
            key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':'+str(min(e[0],d[0])-1)+':'+str(e[1])+':'+str(max(e[0],d[0])-1)+':'+str(e[1]);
            key = supInfo[1]+':'+str(ce[1])+':'+str(min(e[0],d[0])-1)+':'+str(max(e[0],d[0])-1);
            if supInfo[2]=='+': ## positive strand. alt-3 event
              if key in ss3Events: ## duplicate
                dupSS3 += 1;
              else:
                ss3Events[key]=1;
                num3 += 1;
                oFile_3.write(str(num3)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(max(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
            else: ## neg strand. alt-5 event
              if key in ss5Events: ## duplicate
                dupSS5 += 1;
              else:
                ss5Events[key]=1;
                num5 += 1;
                oFile_5.write(str(num5)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(max(e[0],d[0])-1)+'\t'+str(e[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');


      ### Alternative First Exon and Alternative Last Exon

      for tID in genes[gID]: ## examine each transcript to see if the given exon is in a given transcript (index should be 1 or len-2)
        tLen = len(genes[gID][tID]); ## length of a transcript
        if tLen<2: ## not enough exons, skip this
          continue; ## process next transcript

        if [ce[0],ce[1]] == genes[gID][tID][1]: ## current exon is the 2nd exon in the tx
          fEx = genes[gID][tID][0]; ## firstExon, find other tx with different fEX (non overlapping)
          if fullyContainedInInternalExon(fEx, gID): ## fEx is fully contained.. process next transcript
            continue;
          for ctxID in genes[gID]: ## need to examine each tx, candidate transcript id
            if len(genes[gID][ctxID])<2: ## not enough exons, skip this
              continue; ## process next candidate transcript
            if [ce[0],ce[1]] == genes[gID][ctxID][1]: ## the target exon is the 2nd exon in the ctx
              cfEx = genes[gID][ctxID][0]; ## candidate first exon
              if cfEx[0]>fEx[1] or fEx[0]>cfEx[1] : ### non-overlapping exon with smaller coord

                #### should not fully contained in an internal exon of other transcript
                if fullyContainedInInternalExon(cfEx, gID): ## cfEx is fully contained.. process next candidate transcript
                  continue;
  
                if supInfo[2]=='+': ## positive strand. AFE, alt first exon event
                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
                  if key in afeEvents: ## already have this one..
                    pass; ## do nothing
                  else: ## new AFE
                    afeEvents[key] =1;
                    numAFE += 1;
                    oFile_afe.write(str(numAFE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
                else: ## neg strand. ALE, alt last exon event
                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
                  if key in aleEvents: ## already have this one..
                    pass; ## do nothing
                  else: ## new ALE
                    aleEvents[key] =1;
                    numALE += 1;
                    oFile_ale.write(str(numALE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');


        if [ce[0],ce[1]] == genes[gID][tID][-2]: ## current exon is the 2nd to the last exon in the tx
          fEx = genes[gID][tID][-1]; ## lastExon, find other tx with different fEX (non overlapping)
          if fullyContainedInInternalExon(fEx, gID): ## fEx is fully contained.. process next transcript
            continue;
          for ctxID in genes[gID]: ## need to examine each tx, candidate transcript id
            if len(genes[gID][ctxID])<2: ## not enough exons, skip this
              continue; ## process next candidate transcript
            if [ce[0],ce[1]] == genes[gID][ctxID][-2]: ## the target exon is the 2nd exon in the ctx
              cfEx = genes[gID][ctxID][-1]; ## candidate last exon
              if cfEx[0]>fEx[1] or fEx[0]>cfEx[1] : ### non-overlapping exon with smaller coord

                #### should not fully contained in an internal exon of other transcript
                if fullyContainedInInternalExon(cfEx, gID): ## cfEx is fully contained.. process next candidate transcript
                  continue;

                if supInfo[2]=='-': ## negative strand. AFE, alt first exon event
                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
                  if key in afeEvents: ## already have this one..
                    pass; ## do nothing
                  else: ## new AFE
                    afeEvents[key] =1;
                    numAFE += 1;
                    oFile_afe.write(str(numAFE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');
                else: ## pos strand. ALE, alt last exon event
                  key = supInfo[1]+':'+str(ce[0]-1)+':'+str(ce[1])+':';
                  key += str(min(fEx[0],cfEx[0])-1)+':'+str(min(fEx[1],cfEx[1]))+':';
                  key += str(max(fEx[0],cfEx[0])-1)+':'+str(max(fEx[1],cfEx[1]));
                  if key in aleEvents: ## already have this one..
                    pass; ## do nothing
                  else: ## new ALE
                    aleEvents[key] =1;
                    numALE += 1;
                    oFile_ale.write(str(numALE)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(max(fEx[0],cfEx[0])-1)+'\t'+str(max(fEx[1],cfEx[1]))+'\t'+str(min(fEx[0],cfEx[0])-1)+'\t'+str(min(fEx[1],cfEx[1]))+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');


        ### Retained Intron events
        for i in range(0,len(uf)): ### going through the upstream flanking exons
          f1 = uf[i]; ## first flanking exon
          for tID in genes[gID]: ## examine each transcript
            if [f1[0],ce[1]] in genes[gID][tID]: ## there is an exon starts from f1 ends at ce, it is retained intron
              #key=supInfo[1]+':'+str(f1[0]-1)+':'+str(ce[1])+':'+str(f1[0]-1)+':'+str(f1[1])+':'+str(ce[0]-1)+':'+str(ce[1]);
              key=supInfo[1]+':'+str(f1[1])+':'+str(ce[0]-1);
              if key in riEvents: ## already have this skipping events
                dupRI +=1;
                continue; ## next transcript
              else: ## new key, write it
                riEvents[key]=1;
                numRI += 1;
                oFile_ri.write(str(numRI)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(f1[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\n');

        for i in range(0,len(df)): ### going through the downstream flanking exons
          f1 = df[i]; ## first flanking exon
          for tID in genes[gID]: ## examine each transcript
            if [ce[0],f1[1]] in genes[gID][tID]: ## there is an exon starts from ce ends at f1, it is retained intron
              key=supInfo[1]+':'+str(ce[1])+':'+str(f1[0]-1);
              if key in riEvents: ## already have this skipping events
                dupRI +=1;
                continue; ## next transcript
              else: ## new key, write it
                riEvents[key]=1;
                numRI += 1;
                oFile_ri.write(str(numRI)+'\t'+gID+'\t'+'\t'.join(supInfo)+'\t'+str(ce[0]-1)+'\t'+str(f1[1])+'\t'+str(ce[0]-1)+'\t'+str(ce[1])+'\t'+str(f1[0]-1)+'\t'+str(f1[1])+'\n');

  ### end of else: end of merging genes
### end for gID in genes   

logging.debug("Done processing each gene from dictionary to compile AS events");
#
logging.debug("Found %d exon skipping events, from dic %d" % (numSkippingEvents, len(sEvents)));
logging.debug("duplicate skipping events: %d" % dupSE);
#
logging.debug("Found %d exon MX events, from dic %d" % (numMXEvents, len(mxEvents)));
logging.debug("duplicate MXE events: %d" % dupMXE);
#logging.debug("Filtered MXE events: %d" % filteredMXE);
#
logging.debug("Found %d alt SS events" % (len(ss3Events)+len(ss5Events)));
logging.debug("There are %d alt 5 SS events and %d alt 3 SS events." % (num5,num3));
logging.debug("duplicate alt-5 SS events: %d" % dupSS5);
logging.debug("duplicate alt-3 SS events: %d" % dupSS3);
#
#
logging.debug("Found %d AFE events, from dic %d" % (numAFE, len(afeEvents)));
#logging.debug("duplicate AFE events: %d" % dupAFE);
#
logging.debug("Found %d ALE events, from dic %d" % (numALE, len(aleEvents)));
#logging.debug("duplicate ALE events: %d" % dupALE);
#
logging.debug("Found %d RI events, from dic %d" % (numRI, len(riEvents)));
logging.debug("duplicate RI events: %d" % dupRI);
#

iFile.close();
oFile_3.close();
oFile_5.close();
oFile_ce.close();
oFile_mxe.close();
#oFile_mxe_filtered.close();
oFile_afe.close();
oFile_ale.close();
oFile_ri.close();

#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
