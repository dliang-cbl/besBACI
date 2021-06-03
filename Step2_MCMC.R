rm(list=ls())
library(callr)
## MCMC parameters
nchains <- 2
iters <- 10
burnin <- 5
thin <- 1

for(i in 1:nchains){
  for(x in c("Q","TN","TP")){
    res <- paste0("MCMC",x,"B",i,".Rout")
    cmd_ <- paste0('RCMD BATCH "--vanilla --slave --args ',
                   x,' ',i,' ',iters,' ', burnin,' ',thin,
                   '" MCMCInternal.R ',res)
    cat(cmd_,'\n')
    rcmd_bg(cmd=cmd_)
    Sys.sleep(1)
  }
}

