#!/usr/bin/env Rscript --slave --vanilla 

#############
# sQTL signal  of one single exon overlapping GWAS signal


args = commandArgs(TRUE)

if(length(args) < 4){
    cat("No arguments supplied. Usage:\n")
    helpstr = "sQTLGWAS.oneexon.R exoninforfile gwascatalogfile geneannotation.gtf  exonindex [pvaluecutoff]
## pvaluecutoff defines the signficance cutoff of GLiMMPS for plotting, which is optional, default = 1e-5

Example :
sQTLGWAS.oneexon.R alltype/bychrs/exonsinfor.plink.5reads.chr2.txt  ../gwascatalog.txt ~/Research/sQTL/GLiMMPScode/reference/Ensembl_r65.gtf.gff 1117 3.7e-6
\n"
    cat (helpstr)
    quit(status=1)

}

exoninforfile <-  "alltype/bychrs/exonsinfor.plink.5reads.chr2.txt" 
gwasfile <-"../gwascatalog.txt"
gtffile <- "~/Research/sQTL/GLiMMPScode/reference/Ensembl_r65.gtf.gff"
exonindex <-  1117 

pval.cutoff<-  10^(-5) ## default pvalue cutoff for signfiance of GLiMMPS
## pval.cutoff <- 3.7*10^(-6) 

PLINKPATH="~/bin/"

exoninforfile = args[1]
gwasfile = args[2]
gtffile= args[3]
exonindex = as.numeric(args[4])


if(length(args) >=5 ){
pval.cutoff<-   as.numeric(args[5])
}



cat ("Input parameters:\n")
cat(c("exoninforfile = ",exoninforfile,"\ngwasfile = ",gwasfile,"\ngtffile = ",gtffile, "\nexonindex = " ,exonindex,"\n") ) # 
cat (paste("pvalue cutoff =",pval.cutoff,"\n") )




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

TMPGENODIR <- "tmpgeno"
TMPASSODIR <- "tmpasso"
PLINKPATH <- "" # "~/bin/"

source("GLiMMPS_functions.R")


 
tmpgenoprefix <- paste(TMPGENODIR,"/chr",chr,"/",targetexonID,sep="")


############################################################################################################################################################
# 
##############
# read in genotype information
tmpgenofile <- paste(TMPGENODIR,"/chr",chr,"/",targetexonID,".raw",sep="")
tmpmapfile <- paste(TMPGENODIR,"/chr",chr,"/",targetexonID,".map",sep="")
map <- as.matrix(read.csv(tmpmapfile, header=F,sep="\t"))
colnames(map) <- c("Chr","SNPID","cM","Pos")


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
minsnp <- seq(1,length(pvals.sig[,1]))[(  pvals.sig$pvals.glmm <= pval.cutoff & !is.na( pvals.sig$pvals.glmm)  )]  
if(length(minsnp)>0) {  ## if there is signficant sQTL SNP
  cat(paste("There are ",length(minsnp), "significant SNPs passed the pvalue.cutoff=",pval.cutoff,"!\n"))

  if(length(minsnp)>1) {  minsnp <- minsnp[rank(dist2exon2[minsnp] ,ties.method="first")==1] }  ## significant sQTL SNP that are closest to the exon Splice site. 
  
  mostsig.snpID <- as.character(pvals.sig[minsnp,2])
  
 # psi.geno.plot (y,n,mostsig.snpID,genesymbol, targetexonID,exon.coordinate, outputdir)
}


cat(paste("Splicing SNP:",mostsig.snpID,"\n"))



################################################################################################
# 
# substract the gene annotation 

ensemblID <- as.character(exon.infor[exonindex,2])

onegene.gtffile <- paste(TMPASSODIR,"/chr",chr,"/",ensemblID,".gtf",sep="")
subgtfcmd <- paste ("grep ",ensemblID," ", gtffile,">",onegene.gtffile,sep="")
system(subgtfcmd)

gene.gtf <- read.csv(onegene.gtffile , sep="\t",header=F)
#gene.gtf <- read.csv(paste(HOME,"/ITPA.gtf",sep=""), sep="\t",header=F)

exon.gtf <- gene.gtf[ gene.gtf[,3]=="exon",]
CDS.gtf <- gene.gtf[ gene.gtf[,3]=="CDS",]

transcripts.names <- unique(exon.gtf[,9])

annoexon.indexes <- rep(0, length(transcripts.names))

for ( ti in 1: length(transcripts.names)) {

  subexon.gtf <- exon.gtf [ exon.gtf[,9]== transcripts.names[ti] ,]
  subexon.gtf <- subexon.gtf [order(subexon.gtf[,4] ),]
  oneindex <- seq(1, length(subexon.gtf[,1])) [subexon.gtf[,4]== (exon.start+1) & subexon.gtf[,5]==exon.end ]
  if (length(oneindex)>0) {
    annoexon.indexes[ti] <- oneindex
  }
}
# Exon index in the gene annotation
annoexon.geneindex <- max(annoexon.indexes,na.rm=T)

####################################################
#  gwasSNP infor

GWASsnp.infor <- read.csv(gwasfile,header=T,sep="\t")

gwasSNPs<-as.character(GWASsnp.infor$SNPs)

write.table(unique(gwasSNPs),"gwascatalog_SNPids.txt",col.names=F,row.names=F,quote=F)

LDoutfile <- paste(tmpgenoprefix, "_LDgwasSNP_rsq0.8_200kb",sep="")

plinkcmmd <-  paste(PLINKPATH, "plink --file ",tmpgenoprefix  ,"  --r2   --ld-snp-list gwascatalog_SNPids.txt --ld-window-kb 200  --ld-window 99999  --ld-window-r2 0.8 --out ", LDoutfile,sep="")
print(plinkcmmd)
system(plinkcmmd)


gwasSNPsLD <- read.table(paste(LDoutfile,".ld",sep=""),header=T, as.is=T) ## 

################
# read in the pvalue file
pvalfile <- paste(TMPASSODIR,"/chr",chr,"/",targetexonID,".asso",sep="")
pvals.sig <- read.table(pvalfile,header=T,sep="\t",as.is=T)   

#pval.cutoff <- 3.7*10^(-6) 
sQTLSNPs <- pvals.sig$SNPID[(  pvals.sig$pvals.glmm <= pval.cutoff & !is.na( pvals.sig$pvals.glmm)  )]


gwaspeakSNP.ID<- unique(gwasSNPsLD$SNP_A[is.element(gwasSNPsLD$SNP_B, sQTLSNPs)])


subgwasinfor <- GWASsnp.infor[ match(gwaspeakSNP.ID, GWASsnp.infor$SNPs),]
unique(subgwasinfor$Disease.Trait) #3 traits
#[1] Chronic lymphocytic leukemia Multiple sclerosis           Crohn's disease             

gwasTrait <- as.character(subgwasinfor$Disease.Trait)
gwaspeakSNP.pos <- (subgwasinfor$Chr_pos)

LDSNPs.combine <- unique(c( gwaspeakSNP.ID,gwasSNPsLD$SNP_B ))

cat("\nOverlapping GWAS signal: \n")
print(subgwasinfor)

write.table(subgwasinfor,  paste(TMPASSODIR,"/chr",chr,"/",targetexonID,"_OverlappingGWAS.infor.txt",sep=""),row.names=F,col.names=T, quote=F,sep="\t" )
################################################
# GWAS signal overlapping sQTL figure


######## SP140 ###


#exonindex <- seq(1, length(exonIDs))[ exonIDs==targetexonID]
allexon.infor <-exon.infor

chrstr <- allexon.infor[exonindex, 4]
strand<- allexon.infor[exonindex, 5]

chr<- sub("chr","",chrstr)

print(c(exonindex, targetexonID))

genename <-  as.character(allexon.infor$geneSymbol[exonindex])   
genomebrowser.coor <- as.character(allexon.infor$TargetExon[exonindex])  
Exon.start <- allexon.infor[exonindex,6]
Exon.end <- allexon.infor[exonindex,7] 
upstreamExon.start <- allexon.infor[exonindex,10] 
upstreamExon.end <- allexon.infor[exonindex,11]

downstreamExon.start <- allexon.infor[exonindex,12]
downstreamExon.end <- allexon.infor[exonindex,13]
  


#######################################################
#   pvalue plotting.

## read in the snp annotation file

tmpsnpfile <- paste(TMPASSODIR,"/chr",chr,"/",targetexonID,".snpanno",sep="")

tmpsnp.infor <- read.csv(tmpsnpfile, header=T,sep="\t")
snp.poss<- tmpsnp.infor[,4]
tmpSNPIDs<- as.character(tmpsnp.infor[,2])


pch.tmp2  <- rep(19,length(tmpsnp.infor[,1]))
mostsig.snpID  #Most likely causal sQTL SNP:   rs28445040
pch.tmp2[is.element(tmpsnp.infor[,2] ,mostsig.snpID) ] <- 17  ### most closest to SS signficant sQTL SNP is code as "triangle" 
  
pch.tmp2[ is.element(tmpsnp.infor[,2] , gwaspeakSNP.ID)] <- 8  ### GWAS Peak SNP is code as "*" 
pch.tmp2[!is.element(as.character(tmpsnp.infor[,2]),LDSNPs.combine)] <-1  ## those not in LD or overlapping GWAS SNPs are in empty dots, instead of solid dots

tmpSNPIDs<- as.character(tmpsnp.infor[,2])



pch.col <- rep(grey(0.6), length(tmpsnp.infor[,1]))
pch.col[is.element(as.character(tmpsnp.infor[,2]),LDSNPs.combine)] <- 1  ## those in LD or overlapping GWAS SNPs are in black


pch.col[ is.element(tmpsnp.infor[,2] , gwaspeakSNP.ID)] <- "blue"
pch.col[is.element(tmpsnp.infor[,2] ,mostsig.snpID)] <- "red"


#########
pvals <- pvals.sig$pvals.glmm #Chr     SNPID       Pos  glm, glm.quasibin, glmm, lm, beta in 8 columns

#################
## plotting  
FDR0.1 <- pval.cutoff# 3.7*10^(-6)  #significance cutoff


####################
# Plotting only the neighboring 20kb region


pdf ( paste(TMPASSODIR,"/chr",chr, "/GWASoverlap_20kbpeak","_",genename,"_",targetexonID,".pdf",sep="") ,height=8, width=8)

plotregion <- 20000

xlimt <- c(Exon.start-plotregion, Exon.end+plotregion)/(10^6) # tmpsnp.infor[,3], na.rm=T)

#xlimt <- c(min( c(CDS.gtf[,4], CDS.gtf[,5]),na.rm=T), max( c(CDS.gtf[,4], CDS.gtf[,5]),na.rm=T))

ylimt <-c (0, max(-log10(pvals), na.rm=T))

layout(matrix(1:3,3,1), heights = c(rep(4,1),0.35,rep(1,1)))


par(mar=c(4,5,3,1), bg="transparent", cex.lab=2,font.lab=1,las=1,cex.axis=2,cex.main=2 )
main = genename #  paste("GWAS SNP:",gwaspeakSNP.ID, gwasTrait ) #paste("GWASpeak:",gwaspeakSNP.pos,gwaspeakSNP.ID, gwasTrait )

substar <-  is.element(pch.tmp2,c(17,8))
pchsize <-0.95
if (plotregion==200000) {
pchsize <-0.5
}
### plot sQTL pvalue
plot ( snp.poss[!substar]/(10^6), -log10(pvals[!substar]), xlab="Position (Mb)", ylab=expression(paste(-phantom(),"log"[10],"(P)")), main= main, pch= pch.tmp2[!substar] , col=pch.col[!substar], xlim=xlimt, ylim=ylimt, cex=pchsize)
#abline(h= -log10(typeIerror0.01), lty=2)
abline(h= -log10(FDR0.1), lty=2,col=1)
abline(v=Exon.start/(10^6) ,lty=2, col=2,lwd=0.5)
abline(v=Exon.end/(10^6) ,lty=2, col=2,lwd=0.5)
# abline(v=Exon.start-2000 ,lty=2,lwd=0.5)
# abline(v=Exon.end + 2000 ,lty=2, lwd=0.5)

par("new"=T)
plot ( snp.poss[substar]/(10^6), -log10(pvals[substar]), xlab="", ylab="", main= "", pch= pch.tmp2[substar] ,col=pch.col[substar],  xlim=xlimt, ylim=ylimt, cex=pchsize*3,xaxt="n", yaxt="n")

#### add GWAS SNP ID
if (length(gwaspeakSNP.ID)>0) {
text (gwaspeakSNP.pos/(10^6) ,-log10(pvals[is.element(snp.poss,gwaspeakSNP.pos) ]) -0.6 , labels= gwaspeakSNP.ID , adj=c(0,0.5),col=4,cex=2)

legend("bottomleft" , legend= paste(gwaspeakSNP.ID, ":", gwasTrait,sep=" "),  cex= 2, text.col= 4, bg = 'gray90' ,inset = .05 )
}

## Splicing SNP ID
text ( tmpsnp.infor[pch.tmp2==17,4]/(10^6) +diff(xlimt)/100   , -log10(pvals[pch.tmp2==17]) -1.25 , labels= as.character(tmpsnp.infor[pch.tmp2==17,2])  ,xpd=TRUE, adj=c(0,0.5),col=2,cex=2)

legend(xlimt[1], ylimt[2]*0.8 , legend=c("Splicing SNP",paste("GWAS SNP",sep="")), pch=c(17,8) , col = c("red","blue"), cex=1.5, pt.cex= 0.95*3 )


zoomregion<- max(abs(Exon.start-upstreamExon.start),  abs(downstreamExon.end-Exon.end)) +1000
#zoomregion<-3000

####################
# zoom in lines

main= ""# paste(genename, targetexonID) # genomebrowser.coor,
par(mar=c(0,5,0,1))
plot(x=NULL, y=NULL, yaxt="n", xaxt="n", xlab=NA, ylab=NA,   xlim=xlimt,      ylim=c(-1,1), frame.plot =FALSE, main=main, cex=0.5)
lines(c(xlimt[1], (Exon.start-zoomregion)/(10^6)), c(-1,1),lty=2)
lines(c(xlimt[2], (Exon.end+zoomregion)/(10^6)), c(-1,1),lty=2)
# 

# ################
# # plot SNPIDs,
# 
# plot(x=NULL, y=NULL, yaxt="n", xaxt="n", xlab=NA, ylab=NA,   xlim=xlimt,      ylim=c(0,1.2), frame.plot = FALSE, main=main, cex=0.5)
# lines( x=xlimt, y=c(0,0) )
# points (  snp.poss, rep( 0, length(tmpsnp.infor[,1])), pch="|",col=pch.col)
# 
# # text (  tmpsnp.infor[,3], rep( 0.2, length(tmpsnp.infor[,1])), tmpSNPIDs, srt=90,adj=0,cex=0.8, col=pch.col)
# legend("topright", legend= c("HAPMAP","1000G"), text.col=c(2,1),bty="n",cex=1.5)
# 
### plot exon structure ##

plotregion <- zoomregion

xlimt <- c(Exon.start-plotregion, Exon.end+plotregion) # tmpsnp.infor[,3], na.rm=T)
ylimt <-c (-2,2)

main= "" # paste(genename, genomebrowser.coor) # targetexonID) # genomebrowser.coor,
par(mar=c(0,5,0,1))
plot(x=NULL, y=NULL, yaxt="n", xaxt="n", xlab=NA, ylab=NA,   xlim=xlimt,      ylim=c(-1,0.8), frame.plot =FALSE, main=main, cex=0.5)
#lines(c(xlimt[1], Exon.start-2000), c(-1,1),lty=2)
#lines(c(xlimt[2], Exon.end+2000), c(-1,1),lty=2)
#


#################
# plot all exon structure of the gene
exonheight<-0.35 

#if ( xlimt[2] < max( c(CDS.gtf[,4], CDS.gtf[,5])) ) { text ( xlimt[2]   ,0, labels=c("//") ,xpd=TRUE, adj=c(0,0.5)) }

#text ( upstreamExon.start-250   ,0, labels=c("//") ,xpd=TRUE, adj=c(1,0.5),cex=2)
# text ( xlimt[2]   ,0, labels=c("//") ,xpd=TRUE, adj=c(0,0.5),cex=2)

############
# add GWAS SNP label
#lines( x=c(upstreamExon.start-250, downstreamExon.end+250 ), y=c(0,0), col=1,lwd=2)
#lines( x=c(xlimt[1], upstreamExon.start-350 ), y=c(0,0), col=1,lwd=2,lty=2)
#lines( x=c(downstreamExon.end+350, xlimt[2] ), y=c(0,0), col=1,lwd=2,lty=2)


lines( x=c(xlimt[1], upstreamExon.start ), y=c(0,0), col=1,lwd=2)
lines( x=c(upstreamExon.start, xlimt[2] ), y=c(0,0), col=1,lwd=2)
text ( xlimt[1]   ,0, labels=c("//") ,xpd=TRUE, adj=c(1,0.5),cex=2)
text ( xlimt[2]   ,0, labels=c("//") ,xpd=TRUE, adj=c(0,0.5),cex=2)

# text ( xlimt[1] +100   ,0, labels=c("|") ,xpd=TRUE, adj=c(0,0.5), cex=2.5, col=4)
# text ( xlimt[1] +100   ,0.35, labels= gwaspeakSNP.ID[1] ,xpd=TRUE, adj=c(0,0.5),col=4,cex=2)
# text ( xlimt[1] +130   ,-0.05, labels= paste(round( (gwaspeakSNP.pos[2]-gwaspeakSNP.pos[1])/1000,1),"Kb")  ,xpd=TRUE, adj=c(0,1),col=1,cex=2) # "18.5-3 Kb"
# 
# text ( xlimt[1] +1300   ,0, labels=c("|") ,xpd=TRUE, adj=c(0,0.5), cex=2.5, col=4)
# text ( xlimt[1] +900   ,-0.35, labels= gwaspeakSNP.ID[2] ,xpd=TRUE, adj=c(0,0.5),col=4,cex=2)
# text ( xlimt[1] +1350   ,0.05, labels= paste(round( (upstreamExon.start-gwaspeakSNP.pos[2])/1000,1),"Kb")  ,xpd=TRUE, adj=c(0,0),col=1,cex=2) # "3 Kb"
# 
# text ( xlimt[2] -100   ,0, labels=c("|") ,xpd=TRUE, adj=c(0,0.5), cex=2.5, col=4)
# text ( xlimt[2] -130   ,0.35, labels= gwaspeakSNP.ID[3] ,xpd=TRUE, adj=c(0.5,0.5),col=4,cex=2)
# text ( xlimt[2] -120   ,-0.2, labels= paste(round( abs(downstreamExon.end-gwaspeakSNP.pos[3])/1000,1),"Kb")  ,xpd=TRUE, adj=c(1,1),col=1,cex=2) # "18.5 Kb"
# 

### Add splicing SNP label
text ( tmpsnp.infor[pch.tmp2==17,4]    ,exonheight, labels=c("|") ,xpd=TRUE, adj=c(0,0.5), cex=2.5, col=2)
text ( tmpsnp.infor[pch.tmp2==17,4]-250   ,exonheight+0.42, labels= as.character(tmpsnp.infor[pch.tmp2==17,2])  ,xpd=TRUE, adj=c(0,0.5),col=2,cex=2)



#points ( tmpsnp.infor[pch.tmp2==17,4], exonheight, pch="|",col=1, cex=2)

# text (  tmpsnp.infor[,3], rep( 0.2, length(tmpsnp.infor[,1])), tmpSNPIDs, srt=90,adj=0,cex=0.8, col=pch.col)

###############
# plot the target exon and neighboring 2 exons


rect(xleft=Exon.start, xright=Exon.end, ybottom=-exonheight, ytop=exonheight, col=2, border = 2,lwd=1)

rect(xleft=upstreamExon.start, xright=upstreamExon.end, ybottom=-exonheight, ytop=exonheight, col=1, border = 1,lwd=1)

rect(xleft=downstreamExon.start, xright=downstreamExon.end, ybottom=-exonheight, ytop=exonheight, col=1, border = 1,lwd=1)

text ( (Exon.start+Exon.end)/2, -0.8, paste("Exon",annoexon.geneindex), col=2,cex=2) # Exon 7
text ( (upstreamExon.start+upstreamExon.end)/2, -0.8, paste("Exon",annoexon.geneindex-1), col=1,cex=2) # Exon 6
text ( (downstreamExon.start+downstreamExon.end)/2, -0.8, paste("Exon",annoexon.geneindex+1), col=1,cex=2) # Exon 8

### UJ
lines( x=c(upstreamExon.end, (Exon.start+upstreamExon.end)/2  ), y=c(exonheight,0.75), col=1,lwd=1)
lines( x=c(Exon.start, (Exon.start+upstreamExon.end)/2 ) , y=c(exonheight,0.75), col=1,lwd=1)
## DJ
lines( x=c(Exon.end, (Exon.end+downstreamExon.start)/2  ), y=c(exonheight,0.75), col=1,lwd=1)
lines( x=c(downstreamExon.start, (Exon.end+downstreamExon.start)/2 ) , y=c(exonheight,0.75), col=1,lwd=1)
## SJ
lines( x=c(upstreamExon.end, (upstreamExon.end+downstreamExon.start)/2  ), y=c(-exonheight,-0.75), col=1,lwd=1,lty=1)
lines( x=c(downstreamExon.start, (upstreamExon.end+downstreamExon.start)/2 ) , y=c(-exonheight,-0.75), col=1,lwd=1,lty=1)

#lines( x=c(Exon.start-2000, Exon.start), y=c(0,0), col=4)
#lines( x=c(Exon.end, Exon.end+2000), y=c(0,0),col=4)

#rect( xleft=Exon.start-plotregion,xright= Exon.start, ybottom=-0.01,ytop=0.01, col=4,border=4)
#rect( xleft=Exon.end,xright= Exon.end+plotregion, ybottom=-0.01,ytop=0.01, col=4,border= 4)


dev.off()
