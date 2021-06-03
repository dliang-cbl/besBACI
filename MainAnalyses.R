## This code documents the post-MCMC analyses
rm(list=ls())
## load Processed data
load(file="besExampleDump.RData")


## Load MCMC results
library(abind)
source("utils.R")
Q_MCMC <- loadMCMC(list.files(pattern="MCMCQB[0-9]+\\.rda",full.names = T))
TN_MCMC <- loadMCMC(list.files(pattern="MCMCTNB[0-9]+\\.rda",full.names = T))
TP_MCMC <- loadMCMC(list.files(pattern="MCMCTPB[0-9]+\\.rda",full.names = T))

## MCMC Diagnostics - Discharge
plot(ts(Q_MCMC$Tausq[,,1]),
     col=1:3,main="Scale Parameters")

plot(ts(Q_MCMC$Bsub[,,1]),main="Selected Mean Parameters")


fitQ <- fittedMCMC(Q_MCMC)
diagplot(log10(discharge),fitQ,1,ylab="Discharge")

## MCMC Diagnostics - total N
plot(ts(TN_MCMC$Tausq[,,1]),col=1:3,main="Scale Parameters")
plot(ts(TN_MCMC$Bsub[,,1]),col=1:3,main="Selected Mean Parameters")
fitTN <- fittedMCMC(TN_MCMC)
diagplot(log10(totalN),fitTN,1,ylab="totalN")

## MCMC Diagnostics - total P
plot(ts(TP_MCMC$Tausq[,,1]),col=1:3,main="Scale Parameters")
plot(ts(TP_MCMC$Bsub[,,1]),col=1:3,main="Selected Mean Parameters")
fitTP <- fittedMCMC(TP_MCMC)
diagplot(log10(totalP),fitTP,1,ylab="totalP")
