#!/usr/bin/env Rscript --slave --vanilla 

#############
# sQTL analysis for SNPs within 200kb window of one single exon
# Last update: 06/19/2019 (by Yungang Xu, yungang.xu@hotmail.com)

## Arguments needed for program listed at the command line
# exoninforfile
# alljunctionfile
# IJfile
# genoplinkfile
# exonindex
# pval.cutoff [optional]
# maf.cutoff [optional]
source("/mnt/isilon/xing_lab/xuy5/GLIMMPS/Geuvadis/GLiMMPScode/Rscripts/GLiMMPS_functions.new.R")

args = commandArgs(TRUE)
options(warn=-1)
if(length(args) < 5){
    cat("No arguments supplied. Usage:\n")
    helpstr = "sQTLregress.oneexon.R exoninforfile  alljunctionfile InclusionJunctionfile genoplinkfile  exonindex [pvaluecutoff]
## pvaluecutoff defines the signficance cutoff of GLiMMPS for plotting, which is optional, default = 1e-5

Example :
sQTLregress.oneexon.R alltype/bychrs/exonsinfor.plink.5reads.chr2.txt  alltype/bychrs/plink.5reads.allreads.chr2.txt alltype/bychrs/plink.5reads.IJ.chr2.txt  ../Genotype/HAPMAP1000G.CheungCEU41.chr2 1117 1e-5
\n"
    cat (helpstr)
    ##supply default values
    exonindex = 1
    quit(status=1)

}


exoninforfile = args[1]
alljunctionfile = args[2]
IJfile= args[3]
genoplinkfile= args[4]
exonindex = as.numeric(args[5])

pval.cutoff<-  10^(-5) ## default pvalue cutoff for signfiance of GLiMMPS
maf.cutoff <- 0.05
#if(length(args) >=6 ){
#pval.cutoff<-   as.numeric(args[6])
#}
if(length(args) >=6 ){
maf.cutoff <- as.numeric(args[6])
}


if(maf.cutoff == 0) warning("The minor allele frequency (MAF) cutoff was set to 0, which may lead to unconverged fitting and kill the run.")
# ###########
# # testing
# exoninforfile =  "alltype/bychrs/exonsinfor.plink.5reads.chr2.txt" 
# alljunctionfile =  "alltype/bychrs/plink.5reads.allreads.chr2.txt" 
# IJfile =  "alltype/bychrs/plink.5reads.IJ.chr2.txt" 
# genotypeplinkfile =  "../Genotype/HAPMAP1000G.CheungCEU41.chr2" 
# exonindex =  1117 
# 
# # end testing
# ##########


phenofile.IJ = IJfile
phenofile = alljunctionfile
genofile = genoplinkfile
cat ("You'r rnning sQTLregress.oneexon.R:\n")
cat ("Input parameters:\n")
cat(c("exoninforfile = ",exoninforfile,"\nalljunctionfile = ",phenofile,"\nIJfile = ",phenofile.IJ,"\ngenotypeplinkfile = ",genofile, "\nexonindex = " ,exonindex,"\n") ) #  c(exoninforfile, phenofile, phenofile.IJ, genofile,exonindex  )))
cat (paste("pvalue cutoff =",pval.cutoff,"\n") )
cat (paste("MAF cutoff =",maf.cutoff,"\n") )




################################
#  read in the exon information

exon.infor <- read.csv(exoninforfile,header=T,sep="\t")
exonIDs <- as.character(exon.infor[,1])
chrstr <- exon.infor[exonindex, 4]
chr<- sub("chr","",chrstr)
targetexonID <- exonIDs[exonindex]

SNPstartpos <- exon.infor[exonindex,6] -200000
SNPendpos <- exon.infor[exonindex,7] +200000
exonnames <- exonIDs 

#print (c(exoninforfile, genofile,  phenofile, phenofile.IJ))
#print (exonindex)

print(c(exonindex, targetexonID))
##############

TMPGENODIR <- "rerun.tmpgeno"
TMPASSODIR <- "rerun.tmpasso"
PLINKPATH <- "" # "~/bin/"
library(lme4) 


###############################################################################################################################################
# read in phenotype expression levels in plink format
allreads.data <- read.csv( phenofile ,header=T,sep="\t")
nsample <- dim(allreads.data)[1] 
nexons <-dim(allreads.data)[2] -2

allreads.matrix <- allreads.data[,seq(1,nexons)+2 ]
IDs.pheno <- as.character(allreads.data[,2]) 
#phenonames <- colnames(allreads.matrix)

print(nexons)

IJ.data <- read.csv( phenofile.IJ ,header=T,sep="\t")
IJ.matrix <- IJ.data[,seq(1,nexons)+2 ]

### create output directory ##
if (!file.exists( TMPGENODIR )) {system(paste ("mkdir",TMPGENODIR)) }
if (!file.exists( TMPASSODIR )) {system(paste ("mkdir",TMPASSODIR)) }

cmmd <- paste("mkdir ",TMPGENODIR,"/chr",chr,sep="")
if (!file.exists( paste(TMPGENODIR,"/chr",chr,sep="")) ) {system(cmmd)}
  

################################################
#### use plink to extract 1000G genotype around +/- 200kb of each exon

 
   tmpgenoprefix <- paste(TMPGENODIR,"/chr",chr,"/",targetexonID,sep="")

   if (!file.exists(paste(TMPGENODIR,"/chr",chr,"/",targetexonID,".raw",sep="") ) ) {
     plinkcmmd <- paste(PLINKPATH, "plink --tfile ", genofile , " --noweb --recode --chr ",chr, " --from-bp ", SNPstartpos," --to-bp ", SNPendpos,  " --out ",tmpgenoprefix,sep="")
     system(plinkcmmd)
     
     tmpmapfile <- paste(TMPGENODIR,"/chr",chr,"/",targetexonID,".map",sep="")
     if (file.exists(tmpmapfile) ){
 
       plinkcmmd <- paste(PLINKPATH, "plink --tfile ", genofile , " --noweb --recodeA --chr ",chr, " --from-bp ", SNPstartpos," --to-bp ", SNPendpos,  " --out ", tmpgenoprefix,sep="")
       print(plinkcmmd)
       system(plinkcmmd)

       plinkcmmd <- paste(PLINKPATH, "plink --file ",tmpgenoprefix  ," --hardy --out "  , tmpgenoprefix ,sep="")
       print(plinkcmmd)
       system(plinkcmmd)

       awkcmmd <- paste("awk '/ALL|TEST/ {print $1,$2,$3,$4,$5,$6,$7,$8}' ",tmpgenoprefix ,".hwe>",tmpgenoprefix,".hwe2",sep="")
       print(awkcmmd)
       system(awkcmmd) 

       cleancmmd <- paste("rm -rf ",tmpgenoprefix,".log ",tmpgenoprefix,"*nosex ", tmpgenoprefix,".hwe", sep="")
       system(cleancmmd)
      
     }
   }
   
  

############################################################################################################################################################
# do the association for all SNPs near the target exon
#
print("start assocation...")

cmmd <- paste("mkdir ",TMPASSODIR,"/chr",chr,sep="")
if (!file.exists( paste(TMPASSODIR,"/chr",chr,sep="")) ) {system(cmmd)}


  pi <- seq(1,nexons)[exonnames==targetexonID]
  cat ( c("Run associations for Exon", pi,exonnames[pi]) )
##############
# read in genotype information
tmpgenofile <- paste(TMPGENODIR,"/chr",chr,"/",targetexonID,".raw",sep="")
tmpmapfile <- paste(TMPGENODIR,"/chr",chr,"/",targetexonID,".map",sep="")
map <- as.matrix(read.csv(tmpmapfile, header=F,sep="\t"))
colnames(map) <- c("Chr","SNPID","cM","Pos")

genoplink <- read.csv(tmpgenofile, header=T,sep=" ", na.strings=c("NA"))
  
IDs.geno <- as.character(genoplink[,2])
IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]
nsnps <- dim(genoplink)[2] -6
sub.geno <- match(IDs.common , IDs.geno)
geno <- as.matrix(genoplink[sub.geno,seq(1, nsnps)+6])
cat(paste(" sQTL analysis for ", nsnps,"SNPs ...\n"))

### only take those individuals with genotypes, and sort the phenotype to be the same order as in the genotype matrix
sub <- match(IDs.common , IDs.pheno)
#IDs.pheno[sub]
    n<- allreads.matrix[sub,pi]
    y <- IJ.matrix[sub,pi]
    

    pvals.glm <- rep(NA,nsnps)
    pvals.glmquasi <- rep(NA,nsnps)
    pvals.glmm <- rep(NA,nsnps)
    pvals.lm <- rep(NA,nsnps)
    pvals.glmmWald <- rep(NA,nsnps)
    snp.mafs <- rep(NA, nsnps)
    betas <-  rep(NA,nsnps) 

###############
#  start association for psi~ SNP for every SNP. 
nSkipped=0
system.time({
for ( gi in  1: nsnps) {         #701:720){ #
 # cat(paste("SNP ", gi,"\n"))
  ### data association
  SNP=  geno[,gi]
  snp.mafs[gi] <- maf(SNP)
  if( snp.mafs[gi] < maf.cutoff){
    nSkipped=nSkipped+1
    #cat(paste("Warning: SNP", colnames(geno)[gi], "has same genotype across individuals. Skipped (",nSkipped,")!\n"))
    next
  }
  onedata <- list(n=n,y=y, SNP=SNP)

  results.glm <- glm.sQTL ( onedata )
  
  pvals.glm[gi] <- results.glm$pval

  results.quasi <- glmquasi.sQTL ( onedata )
  pvals.glmquasi[gi] <- results.quasi$pval

  ###############
  # GLiMMPS method 
  results.glmm <- glmm.sQTL ( onedata )
  pvals.glmm[gi] <- results.glmm$pval
  betas[gi] <- results.glmm$betas[2]
  ############
  results.lm <- lm.sQTL ( onedata )
  pvals.lm[gi] <- results.lm$pval
  
  results.glmmWald <- try(suppressMessages(glmmWald.sQTL( onedata )), TRUE)
  if(isTRUE(class(results.glmmWald)=="try-error")){
    pvals.glmmWald[gi] <- NA
#    cat(paste("AS event:",exonindex, targetexonID," SNP:", colnames(geno)[gi], "!\n"))
    next
  }else{
    pvals.glmmWald[gi] <- results.glmmWald$pval
  }

 }

cat(paste("Warning: ",nSkipped,"SNPs","have the MAFs less than the cutoff",maf.cutoff,"!\n"))

tmpout <- cbind(map[,c(1,2,4)], formatC(pvals.glm,format="e",digits=3),    formatC(pvals.glmquasi,format="e",digits=3),   formatC(pvals.glmm,format="e",digits=3),   formatC(pvals.lm,format="e",digits=3),   formatC(pvals.glmmWald,format="e",digits=3), round(betas,3), round(snp.mafs,4))

colnames(tmpout) <- c("Chr","SNPID","Pos","pvals.glm","pvals.glmquasi","pvals.glmm","pvals.lm","pvals.glmmWald","Beta", "MAF")
write.table(tmpout  , paste(TMPASSODIR,"/chr",chr,"/",targetexonID,".asso",sep=""),row.names=F,col.names=T, quote=F,sep="\t" )
# pvals.glmquasi,pvals.glmm,pvals.lm, pvals.glmmWald

######################
cat("Finished association!... \nChecking signficant sQTL SNPs and plot one signficant SNP that is closest to the target exon Splice Site\n")



############################################################################################################################################################
# check if there is any SNP pass the pval.cutoff, and plot psi against the significant sQTL SNP that is closest to the target exon Splice Site.
#################


# read in the SNP allele information 
SNPinforfile <- paste(tmpgenoprefix,".hwe2",sep="")
SNP.infor <- read.csv(SNPinforfile,header=T,sep=" ")

genesymbol <- as.character(exon.infor[exonindex,3])
exon.coordinate <- as.character(exon.infor[exonindex,17])


#################
# read in the pvalue file
pvalfile <- paste(TMPASSODIR,"/chr",chr,"/",targetexonID,".asso",sep="")
pvals.sig <- read.table(pvalfile,header=T,sep="\t",as.is=T)   # colClasses=list("'integer","character", "integer", "double","double","double","double","double","double") )

###################
# annotate SNPs
exon.start=  exon.infor[exonindex,6]
exon.end =  exon.infor[exonindex,7]
exon.strand =  as.character(exon.infor[exonindex,5])

snpdist2exon <- cbind((pvals.sig$Pos-exon.start), (pvals.sig$Pos-exon.end)) 

dist2exon <- rep(NA, length(pvals.sig[,1]))
dist2exon <- apply(abs(snpdist2exon),1,min,na.rm=T)
dist2exon[(snpdist2exon[,1] > 0) & (snpdist2exon[,2] <= 0)] <- 0 # inside exon, set distance to 0. 

dist2exon2 <- dist2exon

allsnps.anno <- SNP.annotation( pvals.sig$Pos, exon.start, exon.end,exon.strand)
write.table(cbind(map,allsnps.anno ) , paste(TMPASSODIR,"/chr",chr,"/",targetexonID,".snpanno",sep=""),row.names=F,col.names=T, quote=F,sep="\t" )


dist2exon2[allsnps.anno=="5SS" |allsnps.anno=="3SS" |allsnps.anno=="3SSAG" |allsnps.anno=="5SSGT" ]  <- (-9) # at SS, set to (-9).


## significant sQTL SNP that are closest to the exon Splice site. 

#pval.cutoff<-  10^(-5)

outputdir <-  paste(TMPASSODIR,"/chr",chr,"/",sep="")
suppressWarnings(minsnp <- seq(1,length(pvals.sig[,1]))[(as.numeric(pvals.sig$pvals.glmm) <= pval.cutoff & !is.na( as.numeric(pvals.sig$pvals.glmm)))])
if(length(minsnp)>0) {  ## if there is signficant sQTL SNP
  cat(paste("There are ",length(minsnp), "significant SNPs passed the pvalue.cutoff=",pval.cutoff,"!\n"))

  if(length(minsnp)>1) {  minsnp <- minsnp[rank(dist2exon2[minsnp] ,ties.method="first")==1] }  ## significant sQTL SNP that are closest to the exon Splice site. 
  
  mostsig.snpID <- as.character(pvals.sig[minsnp,2])
  
  psi.geno.plot (y,n,mostsig.snpID,genesymbol, targetexonID,exon.coordinate, outputdir)
  cat("Finished plotting!\n")
}else{
  cat(paste("There is no significant SNPs passed the pvalue.cutoff=",pval.cutoff,"!\n"))
}
})

##  psi.geno.plot (y,n,"rs28445040",genesymbol, targetexonID,exon.coordinate)

