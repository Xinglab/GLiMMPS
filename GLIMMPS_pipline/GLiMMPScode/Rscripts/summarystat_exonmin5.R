#!/usr/bin/env Rscript --slave --vanilla 

############
# read in files 

setwd("./")
args <- commandArgs(T)
type <- "alltype"
PROJECTTITLE <- args[1] ####  Change this to your project title in the config.GLiMMPS.txt file

syscmmd <- paste("cp AScounts/",PROJECTTITLE,".",type,".exoninfor.txt", " AScounts/",PROJECTTITLE,".",type,".exoninforshort.txt"    ,sep="")
cat (syscmmd)
system( syscmmd )

CEU.exoninfor <- read.csv( paste("AScounts/",PROJECTTITLE,".",type,".exoninforshort.txt",sep=""),header=T,sep="\t",row.names=1)
 
 AStype <- sapply( strsplit(rownames(CEU.exoninfor),split="_"), "[[",1)
 TargetExon <- paste(CEU.exoninfor[,4],CEU.exoninfor[,3],":",CEU.exoninfor[,5],"-",CEU.exoninfor[,6],sep="")
 AlternativeExon <- paste(CEU.exoninfor[,4],CEU.exoninfor[,3],":",CEU.exoninfor[,7],"-",CEU.exoninfor[,8],sep="")
 AlternativeExon[ AStype == "A5SS"| AStype=="A3SS" ] <- paste(CEU.exoninfor[,4],CEU.exoninfor[,3],":",CEU.exoninfor[,9],"-",CEU.exoninfor[,10],sep="")[ AStype == "A5SS"| AStype=="A3SS" ]
 ExonID <- rownames(CEU.exoninfor)
 write.table(cbind(ExonID, CEU.exoninfor, AStype, TargetExon, AlternativeExon),  paste("AScounts/",PROJECTTITLE,".",type,".exoninfor.txt",sep=""),col.names=T, row.names=F,sep="\t", quote=F)
 


CEU.exoninfor <- read.csv( paste("AScounts/",PROJECTTITLE,".",type,".exoninfor.txt",sep=""),header=T,sep="\t",row.names=1)



#### filter out those with really short exons (mappable length =0) 
validjunc <- (CEU.exoninfor$IncFormLen >0 & CEU.exoninfor$SkipFormLen >0)
CEU.exoninfor<- CEU.exoninfor[validjunc,]

CEU_IJ <- read.csv( paste("AScounts/",PROJECTTITLE,".",type,".IJ.txt",sep=""),header=T,sep="\t",row.names=1)
CEU_IJ <- CEU_IJ[validjunc,]

nsample <- length(CEU_IJ[1,])
nexons <- length(CEU_IJ[,1])

IJ.matrix <- as.matrix(CEU_IJ)


CEU_SJ <- read.csv(paste("AScounts/",PROJECTTITLE,".",type,".SJ.txt",sep=""),header=T,sep="\t",row.names=1)
CEU_SJ <- CEU_SJ[validjunc,]

SJ.matrix <- as.matrix(CEU_SJ) 


#title <- as.matrix(read.csv("title.txt",header=F,sep="\t"))
CEUIDs <- colnames(CEU_IJ)  

if (!file.exists(type) ) {system(paste("mkdir ",type))}

setwd(type)
##########################
# 

AStype <- sapply( strsplit(rownames(CEU.exoninfor),split="_"), "[[",1)
# for type SE and RI,  effective number of Inclusion Junction reads should be divided by 2, because of 2 possible inlcusion junction (upstream and downstream), and only one possible skipping junction.
IJ.matrix[AStype=="SE",]  <- ceiling((IJ.matrix[AStype=="SE",])/2)
IJ.matrix[AStype=="RI",]  <- ceiling((IJ.matrix[AStype=="RI",])/2)

allreads.matrix <- IJ.matrix+SJ.matrix

newpsi <- IJ.matrix/allreads.matrix




#uniquepsi <- function(x) { return (length(unique(x[!is.na(x)]))) }

#psicounts <- apply(newpsi,1,uniquepsi)
maxpsi <- apply(newpsi,1,max,na.rm=T)
minpsi <- apply(newpsi,1,min,na.rm=T)
psirange <- (maxpsi - minpsi) 
psirange[is.na(psirange) | psirange<=0 | psirange>1] <- 0

rownames(newpsi) <- rownames(CEU_SJ)
colnames(newpsi) <- CEUIDs


medianpsi <- apply(newpsi,1,median,na.rm=T)
meanpsi <- apply(newpsi,1,mean,na.rm=T)
meanallreads <- apply(allreads.matrix,1,mean,na.rm=T)
medianallreads <- apply(allreads.matrix,1,median,na.rm=T)

one.matrix <- (newpsi != medianpsi & (!is.na(newpsi)))
onepsi.all <- apply(one.matrix, 1, sum)

############

pdf(paste("hist_juncreads.pdf",sep=""))
par(mfrow=c(2,2))
zero.matrix <- (allreads.matrix ==0)
zerocounts.all <- apply(zero.matrix, 1, sum)
hist(zerocounts.all, xlab="# Individuals", main="Exons with no junction reads")

zero.matrix <- (allreads.matrix >=10)
zerocounts.all <- apply(zero.matrix, 1, sum)
hist(zerocounts.all, xlab="# Individuals", main="Exons with >=10 junction reads")
zero.matrix <- (allreads.matrix >=20)
zerocounts.all <- apply(zero.matrix, 1, sum)
hist(zerocounts.all, xlab="# Individuals", main="Exons with >=20 junction reads")

zero.matrix <- (allreads.matrix >=50)
zerocounts.all <- apply(zero.matrix, 1, sum)
hist(zerocounts.all, xlab="# Individuals", main="Exons with >=50 junction reads")

dev.off()



#############
# max, median, quantile allreads of each gene.
quantile.matrix <- apply(allreads.matrix, 1, quantile, probs = seq(0, 1, 0.1) )
colnames(quantile.matrix) <- CEU_SJ[,1]


quantile.matrix1 <- quantile.matrix
quantile.matrix1[quantile.matrix1>10] <- 10


pdf (paste("histpercentile_log10readcounts.pdf",sep=""), height=8.5, width=11.5)
par(mfrow=c(4,3))
for ( i in 1:11) {
    hist( log10(quantile.matrix[i,]), main = paste(rownames(quantile.matrix)[i], "percentile per exon") , xlab="log10(Total junction reads)" )
}
dev.off()

pdf (paste("histpercentile_readcounts.pdf",sep=""), height=8.5, width=11.5)
par(mfrow=c(4,3))

for ( i in 1:11) {
    hist( (quantile.matrix1[i,]), main = paste(rownames(quantile.matrix)[i], "percentile per exon") , xlab="Total junction reads" )
}

dev.off()

 length(quantile.matrix[6,quantile.matrix[6,]>=10])
 length(quantile.matrix[6,quantile.matrix[3,]>=10])
 length(quantile.matrix[6,quantile.matrix[2,]>=10])
 length(quantile.matrix[6,quantile.matrix[1,]>=10])


 length(quantile.matrix[6,quantile.matrix[6,]>=20])
 length(quantile.matrix[6,quantile.matrix[3,]>=20])
 length(quantile.matrix[6,quantile.matrix[2,]>=20])
 length(quantile.matrix[6,quantile.matrix[1,]>=20])




############
autosomes <- paste("chr",c(seq(1,22), "X"),sep="")
medianallreads <- quantile.matrix[6,]

sub0<- medianallreads>=5  & psirange>0.1 &  is.element(CEU.exoninfor$chr , autosomes[seq(1,22)]) &  onepsi.all >= 3 ## median of total junction reads >=10  ( (UJ+DJ/2 +SJ >=10 ) ; and at least some psi variation range >0.1## 
print(dim(allreads.matrix[sub0,]))
##   exons with median reads >=5  & psirange>0.1  & with >=3 individuals differ from the medianpsivalue 

### reorder by gene position ###
chrs <-as.numeric(sub("chr","",CEU.exoninfor$chr))
chrs[CEU.exoninfor$chr=="chrX"] <-23
chrs[CEU.exoninfor$chr=="chrY"] <-24
chrs[is.na(chrs)] <- 25 # those not on chr 1-22, X and Y.

poss <- chrs*10^9 + CEU.exoninfor[,5] + CEU.exoninfor[,6] - medianallreads/1000 # order by chr, pos, maxReads isoforms first for the same target exon

targetexon <- paste( AStype, CEU.exoninfor[,3] , CEU.exoninfor[,5] , CEU.exoninfor[,6],sep=":") 


sub <- sub0[order(poss)]  &  ( !duplicated( targetexon[order(poss)] ))    # if target exons have multi AS paths, take only the major isoform AS events.

subexoninfor <- CEU.exoninfor[order(poss),] [sub,]

print(dim(subexoninfor))

geneIDs <- as.character(subexoninfor[,1]) 
start200kb <- subexoninfor[,5] -200000
end200kb <- subexoninfor[,6]+200000

if (!file.exists( "allchrs") ) {system("mkdir allchrs")}
if (!file.exists( "bychrs") ) {system("mkdir bychrs")}


write.table( cbind(subexoninfor) ,"allchrs/subexoninfor.txt",row.names=T,sep="\t",quote=F )
write.table(allreads.matrix[order(poss),][sub,],"allchrs/suballreads.txt",row.names=T,sep="\t",quote=F )
write.table(IJ.matrix[order(poss),][sub,],"allchrs/subIJreads.txt",row.names=T,sep="\t",quote=F )
write.table(SJ.matrix[order(poss),][sub,],"allchrs/subSJreads.txt",row.names=T,sep="\t",quote=F )
write.table(newpsi[order(poss),][sub,],"allchrs/subpsi.txt",row.names=T,sep="\t",quote=F )

write.table( cbind(CEU.exoninfor[,1], meanpsi,meanallreads,medianallreads,psirange)[order(poss),][sub,] ,"allchrs/meanpsi_meanN_medianN_psirange.txt",row.names=T,sep="\t",quote=F )


#######
# output by chrs.
# exonifor;  plink format for all junction reads and inclusionjunciton reads

suballreads.matrix <- allreads.matrix[order(poss),][sub,]
subIJ.matrix <- IJ.matrix[order(poss),][sub,]
subSJ.matrix <- SJ.matrix[order(poss),][sub,]


for ( chr in 1:22) {
  
subsetchr  <- ( is.element(subexoninfor$chr, autosomes[chr]) )
ExonID <- rownames(subexoninfor)
write.table( cbind(ExonID, subexoninfor, start200kb,end200kb)[subsetchr,]  ,  paste("bychrs/exonsinfor.plink.5reads.",autosomes[chr],".txt",sep=""),row.names=F, quote=F,sep="\t" )


output <-  t(rbind(CEUIDs,CEUIDs, suballreads.matrix[subsetchr,] ))
colnames(output)[1] <- "FID"
colnames(output)[2] <- "IID"
 dim(output)

write.table(output , paste("bychrs/plink.5reads.allreads.",autosomes[chr],".txt",sep=""),row.names=F, quote=F,sep="\t" )

output <-  t(rbind(CEUIDs,CEUIDs, subIJ.matrix[subsetchr,] ))
colnames(output)[1] <- "FID"
colnames(output)[2] <- "IID"
 dim(output)

write.table(output , paste("bychrs/plink.5reads.IJ.",autosomes[chr],".txt",sep=""),row.names=F, quote=F,sep="\t" )


output <-  t(rbind(CEUIDs,CEUIDs, subSJ.matrix[subsetchr,] ))
colnames(output)[1] <- "FID"
colnames(output)[2] <- "IID"
 dim(output)

write.table(output , paste("bychrs/plink.5reads.SJ.",autosomes[chr],".txt",sep=""),row.names=F, quote=F,sep="\t" )

}


#################################
