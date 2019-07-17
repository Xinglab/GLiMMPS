##################################################
##
## Used to plot the box plot and sashimi plot for a specific AS-sQTL pair
## 
## This R script is supposed to run in a working dir in "Exon_Inc_Simple/", which is the output foder of all GLIMMPs results.
## Parameters:
## <pop> the population ID: CEU,GBR, FIN, YRI,TSI...or the parent folder where to find the subfolder "Exon_Inc_Simple/", which has the GLiMMPS results.
## <chr> the chromosome ID: chr1:chr22
## <AS_ID> the event ID: SE_1
## <SNP_ID> the SNP ID: rs1000000
##
## USAGE: Rscript plotsQTL.R CEU chr1 SE_96902 rs4649563
###################################################
args = commandArgs(TRUE)
pop = args[1]#"CEU
chr = args[2]#"chr1"; chr2
as = args[3] #"SE_96902";SE_43499
snp = args[4] #"rs4649563";rs28445040
outdir = args[5]

if(is.na(outdir)){
  outdir = paste(pop,".sQTL_plots", sep = "")
  if(!dir.exists(outdir)) dir.create(outdir, showWarnings = F)
} 


#if(!require("staplr",character.only = TRUE)) install.packages("staplr", repos = "https://cloud.r-project.org")
#library("staplr",character.only = TRUE)

tmpasso = paste(pop,"/Exon_Inc_Simple/parallel.tmpasso", sep = "")  ## the outputs folder for GLIMMPs sQTL associations (pvalues)
tmpgeno = paste(pop,"/Exon_Inc_Simple/parallel.tmpgeno", sep = "")
#as = "SE_96902"
#pop = "CEU"
etype = strsplit(as,"_")[[1]][1]

eventfile = paste(pop,"/Exon_Inc_Simple/AScounts/",pop,".", etype,".exoninfor.txt", sep = "")

###############################################################################################################################################
########################
# plot functions for GLiMMPS
#
# plot the psi ~ SNP with dot size proportional to total read counts for each individual
psi.geno.plot <- function(y,n,snpID,genesymbol, exonname,exon.coordinate,outputdir ) {
  
  pvals <- pvals.sig[pvals.sig$SNPID==snpID,]
  
  xlabel <- snpID #paste( , "\tMAF=", pvals$MAF,sep="")
  maintitle <- paste(genesymbol,"\n", exon.coordinate,"\np-value=", pvals$pvals.glmm,sep="")
  subtitle <- paste( "(MAF=", pvals$MAF,")",sep="") #exon.chr,":",exon.start,"-",exon.end,sep="")
  title2 <- paste(exonname,genesymbol,snpID,sep="_")
  
  #print(SNP)
  sNum <- table(SNP)
  aaNum <- 0
  abNum <- 0
  bbNum <- 0
  aaNum <- as.numeric(sNum[names(sNum)==0])
  abNum <- as.numeric(sNum[names(sNum)==1])
  bbNum <- as.numeric(sNum[names(sNum)==2])
  aaNum <- max(aaNum, 0)
  abNum <- max(abNum, 0)
  bbNum <- max(bbNum, 0)
  Alleles <- c(paste(geno.0,"(",aaNum,")",sep="")  ,paste(geno.1,"(",abNum,")",sep=""),  paste(geno.2,"(",bbNum,")",sep=""))  # AA/Aa/aa

  N <- length(n)
  

  ######################
  # use the plot function for reach exonpsi-SNP pair
  ####### y, n, SNP, xlabel, title, Alleles (allele information)
  cat ( paste("Generating Figure: ",outputdir,"/psiplot_",title2,".pdf\n",sep=""))
  plot.pdf <- paste(outputdir,"/Boxplot_",title2,".pdf",sep="")
  pdf( plot.pdf, width=4,height=8) 
  
  # par(oma=c(0,0,0,0),mar=c(2.8,4.5,2.8,0.25), cex.lab=1.5,las=1,cex.axis=1.5) #,bty="l")
  par(oma=c(0.5,0.5,0.5,0),mar=c(6,4.5,4,0.5), cex.lab=1.5,las=1,cex.axis=1.5) #,bty="l"), mfrow=c(1,2)
  
  ##############
  ylim.range <- c(0,1) #range(psi,na.rm=T)
  plot(jitter(SNP,factor=0.5), psi,ylab=expression(paste(psi," (RNA-Seq)")), xlab= xlabel,main= maintitle, sub = subtitle, cex.sub=1, xlim=c(-0.25,2.75), xaxt="n",type="n" ,ylim=ylim.range) #ylim= c(0,1))
  points(jitter(SNP,factor=0.5)  ,psi  , pch= 19, cex= log10(n+1)/1 ,col=1)
  mtext(text=Alleles, side=1, at= c(0,1,2),cex=1.5,line= 1)
  #  mtext(text = xlabel, side=1,at=1,cex=1.5,line=1.6)
  
  
  ## legend of dot size 
  points( rep(2.3,6) , seq(1,6)*0.05+0.5  , pch= 19, cex= log10( c(1,5,10,20,50,100)+1)/1  )
  text(rep(2.55,7) , seq(1,7)*0.05+0.5, c(1,5,10,20,50,100,"# reads"))
  
  par("new"=T) # add boxplot on top
  boxplot(psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.75), ylim=ylim.range,at=sort(unique(SNP[!is.na(SNP)])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0.6) ,outline=FALSE)
  
  
  dev.off()
  return(plot.pdf)
}

exoninforfile = paste(pop,"/Exon_Inc_Simple/alltype/bychrs/exonsinfor.plink.5reads.",chr,".txt", sep = "") 
phenofile.IJ = paste(pop,"/Exon_Inc_Simple/alltype/bychrs/plink.5reads.IJ.",chr,".txt", sep = "")
phenofile = paste(pop,"/Exon_Inc_Simple/alltype/bychrs/plink.5reads.allreads.",chr,".txt", sep = "")
#genofile = ../SNPs_plink/CEU.plink.1
SNPinforfile <- paste(tmpgeno,"/", chr,"/", as,".hwe2",sep="")
tmpmapfile <- paste(tmpgeno,"/",chr,"/",as,".map",sep="")
tmpgenofile <- paste(tmpgeno,"/",chr,"/",as,".raw",sep="")

# read in phenotype expression levels in plink format
allreads.data <- read.csv( phenofile ,header=T,sep="\t")
nsample <- dim(allreads.data)[1] 
nexons <-dim(allreads.data)[2] -2
allreads.matrix <- allreads.data[,seq(1,nexons)+2 ]
IDs.pheno <- as.character(allreads.data[,2]) 


genoplink <- read.csv(tmpgenofile, header=T,sep=" ", na.strings=c("NA"))
rownames(genoplink) <- genoplink[,2]
nsnps <- dim(genoplink)[2] -6
IDs.geno <- as.character(genoplink[,2])
IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]
nsnps <- dim(genoplink)[2] -6
sub.geno <- match(IDs.common , IDs.geno)
geno <- as.matrix(genoplink[sub.geno,seq(1, nsnps)+6])


map <- as.matrix(read.csv(tmpmapfile, header=F,sep="\t"))
SNP.infor <- read.csv(SNPinforfile,header=T,sep=" ")
exon.infor <- read.csv(exoninforfile,header=T,sep="\t")
ABs <- as.matrix(SNP.infor)[match(snp, SNP.infor$SNP),c(5,4)]  # A/a (Major/minor)
#snpID = "rs4649563"
gi <- seq(1,length(map[,2]))[as.character(map[,2])==snp ]# the index of the target SNP
SNP <- geno[,gi] ## the list of genotypes across the smaples, {0,1,2}

IJ.data <- read.csv( phenofile.IJ ,header=T,sep="\t")
IJ.matrix <- IJ.data[,seq(1,nexons)+2 ]

n =  allreads.matrix[,colnames(allreads.matrix)==as]
y = IJ.matrix[, colnames(IJ.matrix)==as]

exon.coordinate <- as.character(exon.infor[exon.infor$ExonID==as,17])
genesymbol <- as.character(exon.infor[exon.infor$ExonID==as,3])

## thre real genotypes, represened by ref. or alt. allele
geno.0 = paste(ABs[1],ABs[1], sep = "")
geno.1 = paste(ABs[1],ABs[2], sep = "")
geno.2 = paste(ABs[2],ABs[2], sep = "")


psi <- y/n
psi[n==0] <- 1e-5

# read in the pvalue file
pvalfile <- paste(tmpasso,"/",chr,"/",as,".1exon.parallel.asso",sep="")
pvals.sig <- read.table(pvalfile,header=T,sep="\t",as.is=T)   # colClasses=list("'integer","character", "integer", "double","double","double","double","double","double") )
## plot the boxplot
boxplot.file <- psi.geno.plot (y,n,snp,genesymbol, as,exon.coordinate, outdir)
## generate the code to plot the sashimi plot
## grouped by the genotypes

####------------------------------------------------------------------
##
## Generate code for sashimi plots, using rmats2sashimi (modified)
##
##
##

event.all <- read.table(eventfile, header = T, sep = "\t")
na.fill <- rep("NA", 7) # used to fill the empty cells to make a rMAT output-formated event file

event.tar <- event.all[event.all[,1]==as,]
colnames(event.tar)[1] <- "ID"
event.tar[,3] <- genesymbol
samples <- names(SNP)
bams <- paste(pop,"/bam_unique/", samples,".unique.bam", sep = "")

half.sample <- length(bams)%/%2

bamlist1 <- paste(bams[1:half.sample] , collapse=",")
psilist1 <- paste(psi[1:half.sample], collapse=",")
half.sample <- half.sample+1
bamlist2 <- paste(bams[half.sample:length(bams)] , collapse=",")
psilist2 <- paste(psi[half.sample:length(bams)], collapse=",")

event.tar <- cbind(event.tar, t(as.matrix(na.fill)), psilist1, psilist2)
event.tmp <- paste(outdir, "/", pop, ".", as, ".event.txt", sep = "")
write.table(event.tar, event.tmp, quote = F, row.names = F,sep = "\t")

groupAA <- paste(which(SNP==0) , collapse=",")
groupAB <- paste(which(SNP==1) , collapse=",")
groupBB <- paste(which(SNP==2) , collapse=",")
group.info <- paste(outdir, "/", pop, ".", as, ".grouping.gf", sep = "")

file.create(group.info, overwrite = T)
if(nchar(groupAA)>0) write.table(paste(geno.0,":",groupAA, sep = ""), group.info, quote = F, row.names = F, col.names = F, append = T)
if(nchar(groupAB)>0) write.table(paste(geno.1,":",groupAB, sep = ""), group.info, quote = F, row.names = F, col.names = F, append = T)
if(nchar(groupBB)>0) write.table(paste(geno.2,":",groupBB, sep = ""), group.info, quote = F, row.names = F, col.names = F, append = T)


plot.cmd <- paste("rmats2sashimiplot --b1",bamlist1,"--b2", bamlist2, "--l1 S1 --l2 S2 -o", outdir, "--exon_s 1 --intron_s 5 -t", etype, "-e", event.tmp, "--group-info", group.info, sep=" ")
system(plot.cmd)
system(paste("rm -r ",outdir,"/Sashimi_index* ", outdir, "/Sashimi_plot", sep = "") )
system(paste("mv ",outdir,"/Sashimiplot_", genesymbol,".pdf ", outdir,"/Sashimiplot_",as,"_", genesymbol,"_",snp, ".pdf ", sep=""))
cat( paste("Plot files are in ", outdir, "/\n", sep = ""))
#print(event.tar)
#print(psilist2)
#staple_pdf(NULL, c( boxplot.file, paste(outdir,"/Sashimiplot_", genesymbol,".pdf", sep="")), paste(outdir,"/",pop,'.plots.',boxplot.file, sep=""), T)
#system(paste("rm -r ",outdir,"/Sashimiplot_* ", boxplot.file, sep = "") )















