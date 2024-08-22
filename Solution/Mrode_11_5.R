

rm(list=ls())

source("LM_KU_functions2018.R")




################
# Mrode 11.5.1 #
################

# info from book
animal <- seq(from=13,to=26,by=1)
sire <- c(NA,NA,13,15,15,14,14,14,1,14,14,14,14,14)
dam <- c(NA,NA,4,2,5,6,9,9,3,8,11,10,7,12)
EDC <- c(558,722,300,73,52,87,64,103,13,125,93,66,75,33)
fatDYD <- c(9,13.4,12.7,15.4,5.9,7.7,10.2,4.8,7.6,8.8,9.8,9.2,11.5,13.3)
SNP <- rbind(c(2,0,1,1,0,0,0,2,1,2),c(1,0,0,0,0,2,0,2,1,0),c(1,1,2,1,1,0,0,2,1,2),c(0,0,2,1,0,1,0,2,2,1),c(0,1,1,2,0,0,0,2,1,2),c(1,1,0,1,0,2,0,2,2,1),c(0,0,1,1,0,2,0,2,2,0),c(0,1,1,0,0,1,0,2,2,0),c(2,0,0,0,0,1,2,2,1,2),c(0,0,0,1,1,2,0,2,0,0),c(0,1,1,0,0,1,0,2,2,1),c(1,0,0,0,1,1,0,2,0,0),c(0,0,0,1,1,2,0,2,1,0),c(1,0,1,1,0,2,0,1,0,0))
X <- rep(1,8)
R <- diag(1/EDC[1:8])
Ident <- diag(1,10)
s2a <- 35.241
s2e <- 245

# population info in tables
pop.info <- data.frame(animal=animal,sire=sire,dam=dam,mean=rep(1,length(animal)),EDC=EDC,fatDYD=fatDYD)

# view pop.info
pop.info

# population info in a table, adding info of animal 1-12
pop.info <- data.frame(animal=c(seq(1:12),animal),sire=c(rep(NA,12),sire),dam=c(rep(NA,12),dam),mean=rep(1,length(animal)+12),EDC=c(rep(NA,12),EDC),fatDYD=c(rep(NA,12),fatDYD))

# view pop.info with information on parents
pop.info

# compute allele frequency for the SNPs
p <- apply(SNP,2,sum)/(2*nrow(SNP))

# alpha depends also on allele frequencies
alpha <- 2*sum(p*(1-p))*s2e/s2a

# center genotypes
Z <- SNP - matrix(2*p,nrow(SNP),ncol(SNP),byrow=TRUE)

# R inverse for unweighted version is a diagonal
Rinv <- diag(8)

# phenotypes are fat yield for the reference animals (13-20)
y <- fatDYD[1:8]

# calculate crossprods for the MME (unweighted)
XRinvX <- crossprod(X,Rinv)%*%X
XRinvZ <- crossprod(X,Rinv)%*%Z[1:8,]
ZRinvX <- t(XRinvZ)
ZRinvZ <- crossprod(Z[1:8,],Rinv)%*%Z[1:8,]
XRinvy <- crossprod(X,Rinv)%*%y
ZRinvy <- crossprod(Z[1:8,],Rinv)%*%y

# assemble matrix for MME (unweighted)
MME.LHS <- rbind(cbind(XRinvX,XRinvZ),cbind(ZRinvX,ZRinvZ + Ident*alpha))
MME.RHS <- rbind(XRinvy,ZRinvy)

# finally, obtain the estimates
MME.unweighted <- as.numeric(solve(MME.LHS)%*%MME.RHS)

# view MME.est
MME.unweighted

# obtain R inverse for weighted version
Rinv <- solve(R)

# view Rinv
Rinv

# phenotypes are fat yield for the reference animals (13-20)
y <- fatDYD[1:8]

# calculate crossprods for the MME (weighted)
XRinvX <- crossprod(X,Rinv)%*%X
XRinvZ <- crossprod(X,Rinv)%*%Z[1:8,]
ZRinvX <- t(XRinvZ)
ZRinvZ <- crossprod(Z[1:8,],Rinv)%*%Z[1:8,]
XRinvy <- crossprod(X,Rinv)%*%y
ZRinvy <- crossprod(Z[1:8,],Rinv)%*%y

# assemble matrix for MME (weighted)
MME.LHS <- rbind(cbind(XRinvX,XRinvZ),cbind(ZRinvX,ZRinvZ + Ident*alpha))
MME.RHS <- rbind(XRinvy,ZRinvy)

# finally, obtain the estimates
MME.weighted <- as.numeric(solve(MME.LHS)%*%MME.RHS)

# view MME.est
MME.weighted

# final view with both unweighted and weighted
effects <- data.frame(unweighted=MME.unweighted,unweighted=MME.weighted)
rownames(effects) <- c("mean",paste("SNP_",seq(1:10),sep=""))

# view table
round(effects,3)

# predict individuals 9-14
pred.unweighted <- as.numeric(MME.unweighted[1] + SNP[9:14,]%*%MME.unweighted[2:11])
pred.weighted <- as.numeric(MME.weighted[1] + SNP[9:14,]%*%MME.weighted[2:11])

par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(fatDYD[9:14],pred.unweighted,cex=1.5,pch=c(0,1,2,3,4,5))
abline(lm(pred.unweighted ~ fatDYD[9:14])$coef,lty=2,col=2)
legend("bottomright",pch=c(0,1,2,3,4,5),legend=paste("Animal",seq(from=21,to=26)))
plot(fatDYD[9:14],pred.weighted,cex=1.5,pch=c(0,1,2,3,4,5))
abline(lm(pred.weighted ~ fatDYD[9:14])$coef,lty=2,col=2)
legend("bottomright",pch=c(0,1,2,3,4,5),legend=paste("Animal",seq(from=21,to=26)))

