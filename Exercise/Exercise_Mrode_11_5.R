 

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

