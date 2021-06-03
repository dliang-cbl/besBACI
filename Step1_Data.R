#Maincode_data
rm(list=ls())
bes <- read.csv("besExample.csv")

## define temporal variables
dateR <- as.Date(bes$Date.r,"%m/%d/%Y")
bes$ayear <- format(dateR,"%Y")
bes$amonth <- format(dateR,"%m")
bes$cmon <- with(bes, (as.numeric(ayear)-1998)*12 + as.numeric(amonth))
names(bes)[5] <- "phosphate"

## define site
bes$siteID <- as.numeric(factor(bes$Site))

## manually define adjacency matrix
## assumed all correlated for now
library(spdep)
nsite <- length(unique(bes$siteID))
nb <- vector("list",3)
nb[[1]] <- as.integer(c(2,3))
nb[[2]] <- as.integer(c(1,3))
nb[[3]] <- as.integer(c(1,2))
class(nb) <- "nb"
Q <- -1*nb2mat(nb,style="B",zero.policy=TRUE)
dd <- apply(Q,2,sum)
diag(Q) <- -1*dd

## alternative could try spatial independence
## Q <- diag(nsite)


## helper functions to reshape long format to wide format
## of daily data by site (kg/day)
## data: data frame name
## xvar: site ID
## tvar: time ID
## zvar: response ID
## value: a Matrix with column time and row sites
## i.e. site ordered within time
prepMat <- function(data,xvar,tvar,zvar)
{
  #browser()
  oo <- order(data[,tvar],data[,xvar])
  data <- data[oo,]
  tt <- unique(data[,tvar])
  xx <- unique(data[,xvar])
  mm <- matrix(NA,length(xx),length(tt))
  row.names(mm) <- xx
  colnames(mm) <- tt
  idx <- match(data[,xvar],xx)
  idt <- match(data[,tvar],tt)
  mm[cbind(idx,idt)] <- data[,zvar]
  mm
}

discharge <- prepMat(bes,"fSite","Date","Flow")
totalN <- prepMat(bes,"fSite","Date","totalN")
totalP <- prepMat(bes,"fSite","Date","totalP")


## Covariates
## currently no -covariates
Xmat <- matrix(1,length(c(discharge)),1)

## Internal data structure set up
## mainly for efficient computing

## this could take up to a few hours and 
## > 500 Mb of storage in the working directory
source("m1utils.R")
Qxi <- RW1Q(ncol(totalN))
Qtheta <- seasonalQ(365,ncol(totalN))
eigenQtheta <- eigen(Qtheta,symm=T)
Ptheta <- eigenQtheta$vector
Dtheta <- eigenQtheta$value
eigenQxi <- eigen(Qxi, symm=T)
Pxi <- eigenQxi$vector
Dxi <- eigenQxi$value
eLst <- list(Qxi=Qxi,Qtheta=Qtheta,Ptheta=Ptheta,
             Dtheta=Dtheta,Pxi=Pxi,Dxi=Dxi)

save.image(file="besExampleDump.RData")

prep <- function(){
  rm(list=ls())
  load("loadest_3.rda")
  loadest.df <- subset(loadest.df,Site != "GFCP")
  loadest.df$SiteId <- as.numeric(as.factor(loadest.df$Site))
  loadest.df$fSite <- factor(
    loadest.df$SiteId,labels=c(
      "Baisman Run","Dead Run","Gwynnbrook","Glyndon",
      "Gwynns Run","Villa Nova","McDonogh","Pond Branch"))
  with(loadest.df,table(fSite,Site))
  
  loadestDR2 <- subset(loadest.df, Site %in% c("DRKR","GFGL","POBR"))
  loadestDR2$fSite <- droplevels(loadestDR2$fSite)
  write.csv(loadestDR2,file="besExample.csv",row.names = F)
}