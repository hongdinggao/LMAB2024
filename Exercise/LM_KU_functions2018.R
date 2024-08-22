tmp <- try(library(MASS),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("MASS")
   library(MASS)
} else library(MASS)
tmp <- try(library(stats),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("stats")
   library(stats)
} else library(stats)
tmp <- try(library(base),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("base")
   library(base)
} else library(base)
tmp <- try(library(car),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("car")
   library(car)
} else library(car)
tmp <- try(library(Matrix),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("Matrix")
   library(Matrix)
} else library(Matrix)
tmp <- try(library(regress),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("regress")
   library(regress)
} else library(regress)
tmp <- try(library(pedigree),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("pedigree")
   library(pedigree)
} else library(pedigree)
tmp <- try(library(reshape2),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("reshape2")
   library(reshape2)
} else library(reshape2)
tmp <- try(library(lattice),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("lattice")
   library(lattice)
} else library(lattice)
tmp <- try(library(gridExtra),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("gridExtra")
   library(gridExtra)
} else library(gridExtra)
tmp <- try(library(mvtnorm),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("mvtnorm")
   library(mvtnorm)
} else library(mvtnorm)
tmp <- try(library(scatterplot3d),silent=TRUE)
if(class(tmp) == "try-error")
{
   install.packages("scatterplot3d")
   library(scatterplot3d)
} else library(scatterplot3d)

# # source("/home/au383944/Dropbox/Work/Postdoc/PopSimTools/Functions.R")
# source("/Users/bia/Dropbox/Work/Postdoc/PopSimTools/Functions.R")
# 
# m <- 100
# q <- 10
# 
# # simulate genotypes
# data <- sim.gen.data(2000,m,qtl.pos=sample(1:m,q),MAF=runif(m,0.05,0.5),n.gen=10,R=NULL,genome.out="Zmat")
# dimnames(data$genome) <- list(seq(1:nrow(data$genome)),paste("snp",seq(1:ncol(data$genome)),sep=""))
# W1 <- model.matrix(seq(1:2000) ~ data$pop.info[,3])
# dimnames(W1)[[2]] <- seq(1:ncol(W1))
# W1[,1] <- W1[,1] - apply(W1[,2:ncol(W1)],1,sum)
# W1 <- rbind(matrix(0,sum(is.na(data$pop.info[,3])),ncol(W1)),W1)
# W2 <- model.matrix(seq(1:2000) ~ data$pop.info[,4])
# dimnames(W2)[[2]] <- seq(1:ncol(W2))
# W2[,1] <- W2[,1] - apply(W2[,2:ncol(W2)],1,sum)
# W2 <- rbind(matrix(0,sum(is.na(data$pop.info[,4])),ncol(W2)),W2)
# data$polygenic <- cbind(W1,W2)
# dimnames(data$polygenic) <- list(seq(1:nrow(data$polygenic)),c(paste("p1_",seq(1:ncol(W1)),sep=""),paste("p2_",seq(1:ncol(W2)),sep="")))
# 
# # simulate phenotypes
# pheno1 <- sim.pheno(data$genome,mu=2,0.5)
# pheno2 <- sim.pheno(data$genome[,data$qtl],mu=3,0.4)
# pheno3 <- sim.pheno(list(data$genome,data$polygenic),mu=8,c(0.5,0.1))
# 
# # aggregate phenotypes to data
# data$y1 <- pheno1$y
# data$y2 <- pheno2$y
# data$y3 <- pheno3$y
# 
# # save data
# save(data,file="./DataSets/dataLM1.Rdata")
# save(pheno1,file="./DataSets/pheno1LM1.Rdata")
# save(pheno2,file="./DataSets/pheno2LM1.Rdata")
# save(pheno3,file="./DataSets/pheno3LM1.Rdata")

########################################################################################
# mk1matrix(A)                                                                         #
########################################################################################
# Assembles one single matrix from a list of matrices (that have equal number of rows) #
########################################################################################
#--------------------------------------------------------------------------------------#
# ARGUMENTS                                                                            #
#--------------------------------------------------------------------------------------#
# A            list of matrices                                                        #
#--------------------------------------------------------------------------------------#
# VALUE                                                                                #
#--------------------------------------------------------------------------------------#
# A.all        matrix                                                                  #
########################################################################################

mk1matrix <- function(A)
{
   A.all <- numeric(0)
   for(i in 1:length(A))
   {
      A.all <- cbind(A.all,A[[i]])
   }
   rm(i)
   return(A.all)
}

###########################################################################
# mkped(ID,parent1,parent2)                                               #
###########################################################################
# Obtains the pedigree (A matrix) from the information about a population #
###########################################################################
#-------------------------------------------------------------------------#
# ARGUMENTS                                                               #
#-------------------------------------------------------------------------#
# ID           individuals ID                                             #
# parent1      ID of parent 1                                             #
# parent2      ID of parent 2                                             #
#-------------------------------------------------------------------------#
# VALUE                                                                   #
#-------------------------------------------------------------------------#
# A            pedigree matrix: matrix                                    #
###########################################################################

mkped <- function(ID,parent1,parent2)
{
   ped <- data.frame(ID,parent1,parent2)
   ped[is.na(ped) == TRUE] <- 0
   makeA(ped,rep(TRUE,nrow(ped)))
   A <- read.table("A.txt")
   names(A) <- c("rows","cols","vals")
   Amat1 <- dcast(A, rows~cols, value.var="vals")
   # Amat1 has the first column with the values in A$rows
   # so I did the following to get it as I want
   Amat1 <- Amat1[,-1]
   # create a transpose to make the full matrix
   Amat2 <- t(Amat1)
   # zero it's diagonals so the values are not summed twice
   diag(Amat2) <- rep(0,nrow(Amat2))
   # change NA's to zero to make the sum work
   Amat1[is.na(Amat1)] <- 0
   Amat2[is.na(Amat2)] <- 0
   Amat <- Amat1 + Amat2
   # clear up
   file.remove("A.txt")
   # return pedigree matrix
   A <- as.matrix(Amat)
   return(A)
}



#=========================================================================
# function to set up the A matrix
#=========================================================================
# ARGUMENTS
#-------------------------------------------------------------------------
# ped is a data frame
# Ped  with 3 columns:
# animalID    sireID    damID
#=========================================================================
makeA <- function(ped) {
    nanimal <- nrow(ped)
    A <- matrix(0, nrow=nanimal, ncol=nanimal)
    for (i in 1:nanimal) {
        for (j in 1:i) {
            if (i == j) {
                if (ped[i, 2] != 0 && ped[i, 3] != 0) {
                    A[i, j] <- 1.0 + 0.5*A[ped[i, 2], ped[i, 3]]
                } else {
                    A[i, j] <- 1.0
                }
            } else {
                if (ped[i, 2] != 0 && ped[i, 3] != 0) {
                    A[j, i] <- 0.5*(A[ped[j, 1], ped[i, 2]] + A[ped[j, 1], ped[i, 3]])
                } else if (ped[i, 3] == 0 && ped[i, 2] != 0) {
                    A[j, i] <- 0.5*A[ped[j, 1], ped[i, 2]]
                } else if (ped[i, 2] == 0 && ped[i, 3] != 0) {
                    A[j, i] <- 0.5*A[ped[j, 1], ped[i, 3]]
                } else {
                    A[j, i] <- 0.0
                }
            }
            A[i,j] <- A[j,i]
        }
    }
    return(A)
}



#=========================================================================
# function to directly set up A inverse matrix
# accounting for inbreeding
#=========================================================================
# ARGUMENTS
#-------------------------------------------------------------------------
# ped is a data frame
# Ped  with 4 columns:
# animalID    sireID    damID   inbreedingCoefficients
#=========================================================================
makeAinv <- function(ped) {
    ff <- c()
    for (i in ped[1]) {
        ff[i] = ped[i, 4]
    }
    nanimal <- nrow(ped)
    Ainv <- matrix(0, nrow=nanimal, ncol=nanimal)
    p <- rep(0,3)
    w <- c(1.0, -0.5, -0.5)
    for (i in 1:nanimal) {
        p[1] <- ped[i, 1]
        p[2] <- ped[i, 2]
        p[3] <- ped[i, 3]
        if (p[2] > 0 && p[3] > 0) {
            w[2] <- -0.5
            w[3] <- -0.5
            d <- 4.0/(2 - ff[p[2]] - ff[p[3]])
        } else if (p[2] > 0) {
            w[2] <- -0.5
            w[3] <- 0.0
            d <- 4.0/(3 - ff[p[2]])
        } else if (p[3] > 0) {
            w[2] <- 0.0
            w[3] <- -0.5
            d <- 4.0/(3 - ff[p[3]])
        } else {
            w[2] <- 0.0
            w[3] <- 0.0
            d <- 1.0
        }
        for (i in 1:3)
            if (p[i] > 0) {
                for (j in 1:3) {
                    if (p[j] > 0) {
                        Ainv[p[i], p[j]] <- Ainv[p[i], p[j]] + w[i]*w[j]*d
                    }
                }
            }
    }
    return(Ainv)
}




###########################################################################
# centerGenotypes(M)                                                      #
###########################################################################
# center all genotypes from a genotype matrix (M)                         #
###########################################################################
#-------------------------------------------------------------------------#
# ARGUMENTS                                                               #
#-------------------------------------------------------------------------#
# M            SNP-genotypes: matrix                                      #
#-------------------------------------------------------------------------#
# VALUE                                                                   #
#-------------------------------------------------------------------------#
# Z            centered genotypes: matrix                                 #
###########################################################################

centerGenotypes <- function(M)
{
   Z <- M - matrix(apply(M,2,mean),nrow(M),ncol(M),byrow=TRUE)
   return(Z)
}

###########################################################################
# alleleFreq(M)                                                           #
###########################################################################
# calculate all allele frequencies from a genotype matrix (M)             #
###########################################################################
#-------------------------------------------------------------------------#
# ARGUMENTS                                                               #
#-------------------------------------------------------------------------#
# M            SNP-genotypes: matrix                                      #
#-------------------------------------------------------------------------#
# VALUE                                                                   #
#-------------------------------------------------------------------------#
# p            allele frequencies: vector                                 #
###########################################################################

alleleFreq <- function(M)
{
   p <- apply(M,2,mean)/2
   return(p)
}

##########################################################################################
# singleSNPeffect(y,Z)                                                                   #
##########################################################################################
# Obtains the SNP effects by fitting one SNP a time to the model                         #
##########################################################################################
#----------------------------------------------------------------------------------------#
# ARGUMENTS                                                                              #
#----------------------------------------------------------------------------------------#
# y            a vector containing the phenotypes                                        #
# Z            a design matrix for the random effects                                    #
#----------------------------------------------------------------------------------------#
# VALUE                                                                                  #
#----------------------------------------------------------------------------------------#
# pvals        p-values of the estimated fixed effects                                   #
##########################################################################################

singleSNPeffect <- function(y,Z)
{
   # fit phenotype y with each SNP separately as a fixed effect
   pvals <- numeric(0)
   for(i in 1:ncol(Z))
   {
      pvals <- c(pvals,summary(lm((y - mean(y)) ~ Z[,i] - 1))$coefficients[,4])
   }
   rm(i)
   
   plot(-log(pvals),xlab="SNP",ylab="-log(p-value)",pch=8)
   abline(h=-log(0.05),lty=2)
   points(which(-log(pvals) > -log(0.05)),-log(pvals[-log(pvals) > -log(0.05)]),pch=8,col=2)
   
   return(pvals)
}

##########################################################################################
# BLUP.general(y,fixed,random,var.random,var.res,alpha)                                  #
##########################################################################################
# Obtain BLUP estimates for fixed and random effects                                     #
##########################################################################################
#----------------------------------------------------------------------------------------#
# ARGUMENTS                                                                              #
#----------------------------------------------------------------------------------------#
# y            a vector containing the phenotypes                                        #
# fixed        a design vector/matrix for the fixed effects                              #
# random       a matrix or a list of matrix for the random effects                       #
# var.random   a matrix or a list of matrix of variance matrices for the fixed effects   #
# var.res      the variance matrix for the residuals                                     #
# alpha        vector/scalar ratio s2e/s2r, where s2r is the variance the random effects #
#----------------------------------------------------------------------------------------#
# VALUE                                                                                  #
#----------------------------------------------------------------------------------------#
# list of estimated/predicted fixed and random effects                                   #
##########################################################################################

BLUP.general <- function(y,fixed,random,var.random,var.res,alpha)
{
   if(class(y) != "numeric")
   {
      stop("y must be numeric")
   }
   if(!(class(fixed) == "numeric" | class(fixed) == "matrix"))
   {
      stop("fixed effects must be either numeric or vector")
   }
   if(class(fixed) == "numeric")
   {
      if(length(y) != length(fixed))
      {
         stop("y and fixed effects must have the same number of observations")
      }
      fixed <- matrix(fixed,length(fixed),1)
   } else if(class(fixed) == "matrix")
      {
         if(length(y) != nrow(fixed))
         {
            stop("y and fixed effects must have the same number of observations")
         }
      }
   if(!(class(random) == "matrix" | class(random) == "list"))
   {
      stop("random effects must be either a matrix or a list of matrices")
   }
   if(class(random) == "matrix")
   {
      if(length(y) != nrow(random))
      {
         stop("y and random effects must have the same number of observations")
      }
      W <- random
   } else if(class(random) == "list")
      {
         if(any(unlist(lapply(random,class)) != "matrix"))
         {
            stop("all random effects must be a matrix")
         }
         W <- mk1matrix(random)
      }
   if(class(random) != class(var.random))
   {
      stop("random effects and variances of random effects must be both of the same class (either both matrices, or both lists of matrices)")
   }
   if(class(var.random) == "matrix")
   {
      if(nrow(var.random) != ncol(var.random))
      {
         stop("variance of random effects must be a square matrix")
      }
      if(nrow(var.random) != ncol(random))
      {
         stop("the variance of the random effects must have the same dimension as the number of variables on the random effects matrix")
      }
      if(length(alpha) != 1)
      {
         stop("alpha must have only one element, if the model has only one random effect")
      }
      S <- var.random*(1/alpha)
      Sinv <- try(solve(S),silent=TRUE)
      class.Sinv <- class(Sinv)
      while(class.Sinv == "try-error")
      {
         S <- S + diag(10e-6,nrow(S))
         Sinv <- try(solve(S),silent=TRUE)
         class.Sinv <- class(Sinv)
      }
      rm(S)
      rm(class.Sinv)
   } else if(class(var.random) == "list")
      {
         if(length(random) != length(var.random))
         {
            stop("the number of random effects must be the same as the number of variance matrices for the random efects")
         }
         if(any(unlist(lapply(var.random,nrow)) != unlist(lapply(var.random,ncol))))
         {
            stop("variances of the random effects must all be a square matrix")
         }
         if(any(unlist(lapply(random,ncol)) != unlist(lapply(var.random,ncol))))
         {
            stop("the variances of the random effects must have the same dimension as the number of variables on the random effects matrices")
         }
         if(length(alpha) != length(random))
         {
            stop("alpha must have the same number of elements as random effects")
         }
         Sinv <- list(0)
         for(i in 1:length(var.random))
         {
            S <- var.random[[i]]*(1/alpha[i])
            Sinv.tmp <- try(solve(S),silent=TRUE)
            class.Sinv <- class(Sinv.tmp)
            while(class.Sinv == "try-error")
            {
               S <- S + diag(10e-6,nrow(S))
               Sinv.tmp <- try(solve(S),silent=TRUE)
               class.Sinv <- class(Sinv.tmp)
            }
            Sinv[[i]] <- Sinv.tmp
            rm(S)
            rm(Sinv.tmp)
            rm(class.Sinv)
         }
         rm(i)
         Sinv <- as.matrix(bdiag(Sinv))
      }
   if(nrow(var.res) != ncol(var.res))
   {
      stop("variance of the residuals must be a square matrix")
   }
   if(nrow(var.res) != length(y))
   {
      stop("number of residuals must be the same number of phenotype observations")
   }
   y <- matrix(y,length(y),1)
   X <- fixed
   R <- var.res
   Rinv <- try(solve(R),silent=TRUE)
   class.Rinv <- class(Rinv)
   while(class.Rinv == "try-error")
   {
      R <- R + diag(10e-6,nrow(R))
      Rinv <- try(solve(R),silent=TRUE)
      class.Rinv <- class(Rinv)
   }
   rm(R)
   rm(class.Rinv)
   XRX <- crossprod(X,Rinv)%*%X
   XRW <- crossprod(X,Rinv)%*%W
   WRWSinv <- crossprod(W,Rinv)%*%W + Sinv
   XRy <- crossprod(X,Rinv)%*%y
   WRy <- crossprod(W,Rinv)%*%y
   LHS <- rbind(cbind(XRX,XRW),cbind(t(XRW),WRWSinv))
   LHS.inv <- ginv(LHS)
   RHS <- rbind(XRy,WRy)
   MME.sol <- as.numeric(LHS.inv%*%RHS)
   hat.fixed <- MME.sol[1:ncol(X)]
   hat.random <- MME.sol[(ncol(X)+1):length(MME.sol)]
   final.sol <- list(hat.fixed)
   if(class(random) == "matrix")
   {
      final.sol[[2]] <- hat.random
      final.sol[[3]] <- LHS.inv[-seq(from=1,to=length(final.sol[[1]])),-seq(from=1,to=length(final.sol[[1]]))]
      names(final.sol) <- c("fixed","random","C.random")
   } else if(class(random) == "list")
      {
         hat.random.save <- hat.random
         n.random.eff <- unlist(lapply(random,ncol))
         for(i in 1:length(random))
         {
            final.sol[[i+1]] <- hat.random[1:n.random.eff[i]]
            hat.random <- hat.random[-seq(from=1,to=n.random.eff[i])]
         }
         rm(i)
         final.sol[[length(final.sol)+1]] <- LHS.inv[-seq(from=1,to=length(final.sol[[1]])),-seq(from=1,to=length(final.sol[[1]]))]
         names(final.sol) <- c("fixed",paste("random",seq(from=1,to=length(random)),sep=""),"C.random")
      }
   return(final.sol)
}

SNP.BLUP <- function(y,X,Z,R=NULL,alpha,reference.pop=NULL)
{
   if(is.null(R))
   {
      R <- diag(1,length(y))
   }
   if(is.null(reference.pop))
   {
      reference.pop <- seq(from=1,to=length(y))
   }
   tmp <- BLUP.general(y,fixed=X,random=Z[reference.pop,],var.random=diag(1,ncol(Z)),var.res=R,alpha=alpha)
   predictedBV <- as.numeric(Z[-reference.pop,]%*%tmp$random)
   tmp <- list(fixed=tmp[[1]],SNPeffects=tmp[[2]],predictedBV=predictedBV)
   return(tmp)
}

GBLUP <- function(y,X,genotypes,R=NULL,alpha,GD=FALSE,reference.pop=NULL)
{
   if(is.null(R))
   {
      R <- diag(1,length(y))
   }
   if(is.null(reference.pop))
   {
      reference.pop <- seq(from=1,to=length(y))
   }
   Z <- centerGenotypes(genotypes)
   p <- alleleFreq(genotypes)
   if(GD)
   {
      D <- diag(1/(ncol(genotypes)*2*p*(1-p)))
      G <- Z%*%D%*%t(Z)
   } else
      {
         s.var.Z <- 2*sum(p*(1-p))
         G <- tcrossprod(Z)/s.var.Z
      }
   W <- cbind(diag(1,length(reference.pop)),matrix(0,length(reference.pop),nrow(Z)-length(reference.pop)))
   tmp <- BLUP.general(y,fixed=X,random=W,var.random=G,var.res=R,alpha)
   tmp <- list(tmp[[1]],tmp[[2]],tmp[[2]][-reference.pop])
   names(tmp) <- c("fixed","BV","predictedBV")
   return(tmp)
}

SNP.BLUP.polygenic <- function(y,X,Z,A,R=NULL,alpha,reference.pop=NULL)
{
   if(is.null(R))
   {
      R <- diag(1,length(y))
   }
   if(is.null(reference.pop))
   {
      reference.pop <- seq(from=1,to=length(y))
   }
   W <- cbind(diag(1,length(reference.pop)),matrix(0,length(reference.pop),nrow(Z)-length(reference.pop)))
   tmp <- BLUP.general(y,fixed=X,random=list(Z[reference.pop,],W),var.random=list(diag(1,ncol(Z)),A),var.res=R,alpha=alpha)
   predictedBV <- as.numeric(Z[-reference.pop,]%*%tmp[[2]]) + tmp[[3]][-reference.pop]
   tmp <- list(fixed=tmp[[1]],SNPeffects=tmp[[2]],BV=tmp[[3]],predictedBV=predictedBV)
   return(tmp)
}

#=======================================================================
# R function to compute inbreeding coefficients
# ped is a data frame
#=======================================================================
makeF <- function(ped) {
    nanimal <- nrow(ped)
    F <- numeric(nanimal)
    L <- numeric(nanimal)
    D <- numeric(nanimal)
    F0 <- -1.0
    po <- rep(0, nanimal)
    for (i in 1:nanimal) {
        is <- ped[i,2]
        id <- ped[i,3]
        ped[i,2] <- max(is, id)
        ped[i,3] <- min(is, id)
        if (is == 0 && id != 0) {
            D[i] <- 0.5-0.25*(F0 + F[id])
        } else if (is != 0 && id == 0) {
            D[i] <- 0.5-0.25*(F[is] + F0)
        } else if (is == 0 && id == 0) {
            D[i] <- 0.5-0.25*(F0 + F0)
        } else {
            D[i] <- 0.5-0.25*(F[is] + F[id])
        }
        if (i == 1) {
            if (is == 0 || id == 0) {
                F[i] <- 0.0
            }
        }
        if (i > 1) {
            if (is == 0 || id == 0) {
                F[i] <- 0.0
            } else if (ped[i-1,2] == ped[i,2] && ped[i-1,3] == ped[i,3]) {
                F[i] <- F[i-1]
            } else {
                fi <- -1.0
                L[i] <- 1.0
                j <- i
                while (j != 0) {
                    k <- j
                    r <- 0.5*L[k]
                    ks <- ped[k,2]
                    kd <- ped[k,3]
                    if (ks > 0) {
                        while (po[k] > ks) {
                            k <- po[k]
                        }
                        L[ks] <- L[ks] + r
                        if (ks != po[k]) {
                            po[ks] <- po[k]
                            po[k] <- ks
                        }
                        if (kd > 0) {
                            while (po[k] > kd) {
                                k <- po[k]
                            }
                            L[kd] <- L[kd] + r
                            if (kd != po[k]) {
                                po[kd] <- po[k]
                                po[k] <- kd
                            }
                        }
                    }
                    fi <- fi+L[j]*L[j]*D[j]
                    L[j] <- 0.0
                    k <- j
                    j <- po[j]
                    po[k] <- 0
                }
                F[i] <- fi
            }
        }
    }
    return(F)
}
