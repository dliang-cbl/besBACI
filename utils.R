dointerpolate <- function(v)
{ 
  ## a time serie vector v that contains missing values
  ## output: missing data replaced with interpolated values
  x <- 1: length(v)
  id <- which(is.na(v))
  if( length(id) > 0) {
    tmp <- approx( x[-id], v[-id], xout= id,rule=2)$y
    v[id] <- tmp
  }
  v
}
getinitsy <- function(y,dd=NULL)
{
  ## y a matrix of data with missing values ordered site within time
  ## value the interpolated missing values
  if(!is.null(dd)) dim(y) <- dd
  tmp <- t(apply(y,1,dointerpolate))
  c(tmp)[is.na(y)]
}

loadMCMC <- function(x){
  #x <- list.files(pattern="MCMCQB[0-9]+\\.rda",full.names = T)
  storage <- vector("list",length(x))
  for(i in 1:length(x)){
    cat("loading ",x[i],".\n")
    load(file=x[i])
    storage[[i]] <- r
  }
  Tausq <- lapply(storage,function(x) x$Tausq)
  Tausq <- abind(Tausq, along=3)
  B <- lapply(storage,function(x) x$B)
  B <- abind(B, along=3)
  Bsub <- lapply(storage,function(x) x$Bsub)
  Bsub <- abind(Bsub,along=3)
  DIC <- do.call(rbind,lapply(
    storage, function(x) as.data.frame(x$DIC)))
  Y <- lapply(storage,function(x) x$Y)
  Y <- abind(Y,along=3)
  Ypred <- lapply(storage,function(x) x$Ypred)
  Ypred <- abind(Ypred,along=3)
  elapsed <- do.call(rbind,lapply(
    storage, function(x) x$elapsed
  ))
  Ball <- lapply(storage, function(x) x$Ball)
  Ball <- abind(Ball,along=3)
  list(Tausq=Tausq, B=B, Bsub=Bsub,DIC=DIC,Y=Y,Ypred=Ypred,
       Ball=Ball,elapsed=elapsed)
}

fittedMCMC <- function(r,dd=dim(totalN))
{
  yhat <- apply(r$Y[,1,],1,mean) 
  rr <- yhat
  dim(rr) <- dd
  rr
}

diagplot <- function(y,fit,i,main=NULL,ylab="total N Load (ton/day)")
{
  if(is.null(main)){
    main=row.names(y)[i]
  }
  mm <- cbind(y[i,],fit[i,])
  op <- par(mfrow=c(2,2),mar=c(4,4,4,1))
  plot(ts(mm,start=c(1998,10,15),frequency=365.25),
       plot.type="single",lty=c(1,2),col=c(1,2),ylab=ylab,main=main)
  plot(fit[i,],y[i,],xlab="Fitted",ylab="Observed",
       xlim=range(fit[i,],y[i,],na.rm = T),
       ylim=range(fit[i,],y[i,],na.rm = T))
  abline(0,1,col=2,lwd=2)
  rr <- y[i,]-fit[i,]
  plot(fit[i,],rr,xlab="Fitted",ylab="Residual")
  abline(h=c(-2,0,2),lwd=2)
  qqnorm(rr,main="")
  qqline(rr)
  par(op)
}
