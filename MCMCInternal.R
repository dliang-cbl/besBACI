## Main_2.R
## BES analysis code
rm(list=ls())

## receive input arguments
args__ <- commandArgs(trailingOnly = T)

parm <- args__[1] ## parameter to simulate
batch <- as.integer(args__[2]) ## chain id
iters <- as.integer(args__[3]) ## number of iterations
burnin <- as.integer(args__[4]) ## number of burn in
thin <- as.integer(args__[5])  ## number of thin

stopifnot(args__[1] %in% c("Q","TN","TP"))

cat("processing ",parm,", batch=",batch,".\n")

## load R scripts
source("m1utils.R")
source("m3utils.R")
source("utils_1.R")
fit <- source("fitm3d1z.R")$value
 
## load Processed data
load(file="besExampleDump.RData")

## select data to analyze
## and apply transformation
if(parm=="Q"){
  y <- log10(discharge)
}
if(parm=="TN"){
  y <- log10(totalN)
}
if(parm=="TP"){
  y <- log10(totalP)
}


## MCMC parameters
parmsa <- c(1 ,1 ,1 ,1 ,1)
parmsb <- c(1 ,1, 1 ,1 ,1)
nchains <- 1

## set random seeds
set.seed(batch+Sys.time())

## initial values
initsy <-  getinitsy(y)
q <- runif(nchains)
inits <- numeric()
for( k in 1:nchains) inits <- cbind( inits, qgamma(q[k],parmsa,parmsb))

## Run MCMC
r <- fit(y=c(y),X=Xmat,Qs=Q,nt=ncol(y),per=365,
           parmsa=parmsa,parmsb=parmsb,iters=iters,inits=inits,
           initsy=initsy,nchains=nchains,burnin=burnin,thin=thin,
           storeY=TRUE,eLst=eLst)

fn <- paste0("MCMC",parm,"B",batch,".rda")
save(r,file=fn)

