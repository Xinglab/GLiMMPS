########################
# functions for GLiMMPS
#
# Last update: 06/19/2019 (by Yungang Xu, yungang.xu@hotmail.com)
#




library(lme4) 

################################################################################################################################################
# functions for sQTL associations
########
#
# The sQTL functions have input data as a list(y, n, SNP), and output the p-values and beta estimate for fixed effect (intercept and then the SNP effect)  
##################
# generalized linear model 

glm.sQTL<- function (data) {
  options(warn=2) #

  snp.pval <- 1
  betas <-  c(0,0)
  #psi.fitted 
  if ( length(unique(data$SNP[!is.na(data$y) & data$n>0 & !is.na(data$SNP)] )) >1) {
  
  response <- cbind(data$y, data$n-data$y)
  SNP <- data$SNP
  testglm = try(glm(response ~ SNP, family = binomial),FALSE) #
  #summary(results)
 # summary(results)$coefficients
  #summary(results)$fitted.values
  if (!( inherits(testglm,"try-error") || !testglm$converged)) { # only when there is no complete separation problem and it converges in logistic regression
    results <- testglm
    betas <- summary(results)$coefficients[,1]
    snp.pval <- summary(results)$coefficients[2,4]
  }
  #psi.fitted <- results$fitted.value
}
  options(error = NULL, warn = 0)

  return ( list(betas =betas,  pval = snp.pval ) )  
}


##################
# generalized linear model with quasibinomial family to handle overdispersion; it's performing well. 

glmquasi.sQTL<- function (data) {
  options(warn=2) #

  snp.pval <- 1
  betas <-  c(0,0)
  #psi.fitted 
  if ( length(unique(data$SNP[!is.na(data$y) & data$n>0 & !is.na(data$SNP)] )) >1) {
  nomissing <- (!is.na(data$y) & data$n>0 & !is.na(data$SNP))
  response <- cbind(data$y[nomissing], data$n[nomissing]-data$y[nomissing])
  SNP <- data$SNP[nomissing]
  testglm = try(glm(response ~ SNP, family = quasibinomial),FALSE) #quasibinomial(Overdispersed Binomial Regression).
  #summary(results)
 # summary(results)$coefficients
  #summary(results)$fitted.values
  if (!( inherits(testglm,"try-error") || !testglm$converged)) { # only when there is no complete separation problem and it converges in logistic regression
    results <- testglm
    betas <- summary(results)$coefficients[,1]
    snp.pval <- summary(results)$coefficients[2,4]
  }
  #psi.fitted <- results$fitted.value
}
  options(error = NULL, warn = 0)

  return ( list(betas =betas,  pval = snp.pval ) )  
}

##################
# generalized linear mixed model regression  with Wald's test for fixed effect 
glmmWald.sQTL <- function(data) {
    #######  glmm with Wald's test for fixed effect
    # library(lme4)
  options(warn=2) #

  snp.pval <- 1
  betas <-  c(0,0)

  if ( length(unique(data$SNP[!is.na(data$y) & data$n>0 & !is.na(data$SNP)] )) >1) {
  nomissing <- (!is.na(data$y) & data$n>0 & !is.na(data$SNP))
  response <- cbind(data$y[nomissing], data$n[nomissing]-data$y[nomissing])
  SNP <- data$SNP[nomissing]
  obs <- seq(1,length(SNP))

  testglm = try(  suppressMessages( glmer(response ~ SNP +(1|obs), family=binomial)),silent=TRUE) #invidual-level random effect for overdispersion.

  if ( !( inherits(testglm,"try-error")  )) { # only when it converges
    betas <-fixef(testglm)
    if (!is.na(fixef(testglm)[2])) {
      snp.pval <-  summary(testglm)$coefficients[2,4]
    }
  }
  #psi.fitted <- results$fitted.value
}
  options(error = NULL, warn = 0)

  return ( list(betas =betas,  pval = snp.pval ) )  
}


###############
#  GLiMMPS method 
##################
# generalized linear mixed model regression  with LRT test for fixed effect 
glmm.sQTL <- function(data) {
    #######  glmm 
    # library(lme4)
  #options(warn=2) #this turns all warnings into errors, KZ change

  snp.pval <- 1
  betas <-  c(0,0)

  if ( length(unique(data$SNP[!is.na(data$y) & data$n>0 & !is.na(data$SNP)] )) >1) {
  nomissing <- (!is.na(data$y) & data$n>0 & !is.na(data$SNP))
  response <- cbind(data$y[nomissing], data$n[nomissing]-data$y[nomissing])
  SNP <- data$SNP[nomissing]
   obs <- seq(1,length(SNP))
 
  testglm = try( suppressMessages( glmer(response ~ SNP +(1|obs), family=binomial)),silent=TRUE)  #invidual-level random effect for overdispersion. 
  testglm0 = try(suppressMessages( glmer(response ~ 1 +(1|obs), family=binomial)),silent=TRUE) # 

#  if (!( inherits(testglm,"try-error") || inherits(testglm0,"try-error") )) { # only when it converges ## KZchange
    betas <- try(fixef(testglm), silent=TRUE)
    if (! inherits(betas,"try-error") & !is.na(betas[2])) {
          snp.pval <-  anova(testglm, testglm0)$"Pr(>Chisq)"[2]
        } else betas=c(0, 0)

#  } ## KZchange
  #psi.fitted <- results$fitted.value
} 
#  options(error = NULL, warn = 0) ## KZchange

  return ( list(betas =betas,  pval = snp.pval ) )  
}


#######################################
# linear model 
### use psi =y/n and do linear regression ##
lm.sQTL<- function (data) {
  options(warn=2) #

  snp.pval <- 1
  betas <-  c(0,0)
  #psi.fitted 
  if ( length(unique(data$SNP[!is.na(data$y) & data$n>0 & !is.na(data$SNP)] )) >1) {
  
  psi <- data$y/data$n
  psi[data$n==0] <- NA
  SNP <- data$SNP
  testlm = try(lm(psi ~ SNP),FALSE) #quasibinomial(Overdispersed Binomial Regression).
  #summary(results)
 # summary(results)$coefficients
  #summary(results)$fitted.values
  if (!( inherits(testlm,"try-error") )) { # only when there is no complete separation problem and it converges in logistic regression
    results <- testlm
    betas <- summary(results)$coefficients[,1]
    snp.pval <- summary(results)$coefficients[2,4]
  }
  #psi.fitted <- results$fitted.value
}
  options(error = NULL, warn = 0)

  return ( list(betas =betas,  pval = snp.pval ) )  
}



################################################################################################################################################
SNP.annotation <- function(snppos.vector, exon.start, exon.end,exon.strand) {

  SNPtype <- rep("",length(snppos.vector))
  for ( i in 1:length(snppos.vector) ){
    snppos <- snppos.vector[i]
   ############### #########
  dist.3ss <- ( snppos - exon.start)
  dist.5ss <- ( snppos - exon.end)
  
  type <- ">300bp"
  if (    dist.3ss > (-300) && dist.5ss <= 300 ) {
    type <- "<=300bp"
  }
  
# 5'SS :  9 bases long. [3 bases in exon][6 bases in intron] 
# 3'SS : 23 bases long. [20 bases in the intron][3 base in the exon] 

  

    if (    dist.3ss > 0 && dist.5ss <= 0 ) {
      type <- "exon"
    }
  
    if (  dist.3ss > (-20) && dist.3ss<=3 && exon.strand=="+" )  {
      type <- "3SS"
        if (  dist.3ss ==0 || dist.3ss== -1 )  {
          type <- "3SSAG"
        }
    }
   if (  dist.5ss <= (20) && dist.5ss >(-3) && exon.strand=="-" )  {
     type <- "3SS"
     if (  dist.5ss ==1 || dist.5ss== 2 )  {
       type <- "3SSAG"
     }

    }
  
    if (  dist.5ss > (-3) && dist.5ss<=6 && exon.strand=="+")  {
      type <- "5SS"
      if (  dist.5ss ==1 || dist.5ss== 2 )  {
        type <- "5SSGT"
      }
    }
  if (  dist.3ss > (-6) && dist.3ss<=3 && exon.strand=="-")  {
      type <- "5SS"
      if (  dist.3ss ==0 || dist.3ss== -1 )  {
        type <- "5SSGT"
      }
    }
    SNPtype[i] <- type
  } ## end of each SNP

  return(SNPtype)
}


###############################################################################################################################################
########################
# plot functions for GLiMMPS
#
# plot the psi ~ SNP with dot size proportional to total read counts for each individual
psi.geno.plot <- function(y,n,snpID,genesymbol, exonname,exon.coordinate,outputdir ) {

  ABs <- as.matrix(SNP.infor)[match(snpID, SNP.infor$SNP),c(5,4)]  # A/a (Major/minor)
  Alleles <- c(paste(ABs[1],ABs[1],sep="")  ,paste(ABs[1],ABs[2],sep=""),  paste(ABs[2],ABs[2],sep=""))  # AA/Aa/aa

  xlabel <- snpID
  title <- paste( genesymbol,"\n", exon.coordinate,sep="") #exon.chr,":",exon.start,"-",exon.end,sep="")
  title2 <- paste(genesymbol,exonname,snpID,sep="_")

  gi <- seq(1,length(map[,2]))[as.character(map[,2])==snpID ]
  SNP <- geno[,gi]

  psi <- y/n
  psi[n==0] <- NA
  N <- length(n)

  
  ######################
  # use the plot function for reach exonpsi-SNP pair
  ####### y, n, SNP, xlabel, title, Alleles (allele information)
 cat ( paste("Generating Figure: psiplot_",title2,".pdf\n",sep=""))
  pdf( paste(outputdir,"/psiplot_",title2,".pdf",sep=""),width=4,height=8) 
    
# par(oma=c(0,0,0,0),mar=c(2.8,4.5,2.8,0.25), cex.lab=1.5,las=1,cex.axis=1.5) #,bty="l")
  par(oma=c(0.5,0.5,0.5,0),mar=c(4,4.5,4,0.5), cex.lab=1.5,las=1,cex.axis=1.5) #,bty="l")

##############
  ylim.range <- c(0,1) #range(psi,na.rm=T)
 plot(jitter(SNP,factor=0.5), psi,ylab=expression(paste(psi," (RNA-Seq)")), xlab= xlabel, main=title,xlim=c(-0.25,2.75), xaxt="n",type="n" ,ylim=ylim.range) #ylim= c(0,1))
 points(jitter(SNP,factor=0.5)  ,psi  , pch= 19, cex= log10(n+1)/1 ,col=1)
   mtext(text=Alleles, side=1, at= c(0,1,2),cex=1.5,line= 1)
  #  mtext(text = xlabel, side=1,at=1,cex=1.5,line=1.6)


## legend of dot size 
points( rep(2.3,6) , seq(1,6)*0.05+0.5  , pch= 19, cex= log10( c(1,5,10,20,50,100)+1)/1  )
text(rep(2.55,7) , seq(1,7)*0.05+0.5, c(1,5,10,20,50,100,"# reads"))

   par("new"=T) # add boxplot on top
 boxplot(psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.75), ylim=ylim.range,at=sort(unique(SNP[!is.na(SNP)])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0.6) ,outline=FALSE)


dev.off()

}


########################
## Add a new function for MAF calculation based on the SNP genotypes. (by Yungang Xu, 0618/2019)
## Calculate the allele frequeny from a given list of genotypes, e.g. c(0,0,0,1,2,2,1,1,1,1).
## Arguments
## x  a vector of genotypes, across a group of individuals.It may have only one type of genotype
## af a logical value indicating whether return the allele frequencies for all alleles. By default (af = FALSE), only returns the minor allele requency.
## Value
## a numeric mubmer or vectector contains the allele requencies.

maf <- function(x, af = FALSE){ 
  stat <- as.matrix(table(x))
  N <- sum(stat)*2
  if ("0" %in% row.names(stat)){
    aa <- stat["0",]
  }else{
    aa <-0
  }
  if ("1" %in% row.names(stat)){
    ab <- stat["1",]
  }else{
    ab <-0
  }
  if ("2" %in% row.names(stat)){
    bb <- stat["2",]
  }else{
    bb <-0
  }
  
  freq.a <- (aa*2+ab)/N
  freq.b <- (bb*2+ab)/N
  names(freq.a) <- "a"
  names(freq.b) <- "b"
  if(af){
    return(c(freq.a, freq.b))
  }else{
    return(min(freq.a, freq.b))
  }
  
}


########################

