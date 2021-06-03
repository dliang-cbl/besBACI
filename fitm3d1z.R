## version z:
## include a spatial adjacency matrix instead of nr x nc grid
## and temporal eigen decompositions
function( y , X, Qs, nt, per, parmsa, parmsb, iters, inits, initsy, nchains, burnin , 
          thin, mask = NULL,storeB=F,storeY=F,eLst=NULL)
{
# This function assumes that the data y are collected on a regular lattice
#         at evenly spaced time points
# The data must be ordered by site within time.
# y       -- observed data with missing values
# X       -- observed design matrix 
# Qs      -- spatial adjacency matrix
# nt      -- number of time points
# per     -- periodicity of time series
# parmsa  -- alpha parameters of indep. gamma priors on precisions
# parmsb  -- beta parameters of indep. gamma prior on precisions
# inits   -- initial values for tau
# initsy  -- initial values fof missed data  
# nchains -- number of chains
# burnin  -- number of burnins
# thin   -- store every 'slice' iterations
# mask    -- water mask pixels
# storeB  -- whether raw random effects are stored
# storeY  -- whether predictive distribution are stored
# eLst    -- structure matrices for temporal priors
# ns will be calculated as dimension of Qs
# n will be calculated as ns * nt

timestamp1 <- proc.time()

id <- which(is.na(y))

p <- ncol(X)

if(is.null(eLst)){
  Qtheta <- seasonalQ(per,nt)
  Qxi <- RW2Q(nt)
  eigenQtheta <- eigen(Qtheta,symm=T)
  Ptheta <- eigenQtheta$vector
  Dtheta <- eigenQtheta$value
  eigenQxi <- eigen(Qxi, symm=T)
  Pxi <- eigenQxi$vector
  Dxi <- eigenQxi$value
}
else{
  Qtheta <- eLst$Qtheta
  Qxi <- eLst$Qxi
  Ptheta <- eLst$Ptheta
  Dtheta <- eLst$Dtheta
  Pxi <- eLst$Pxi
  Dxi <- eLst$Dxi
}
#Qphi <- CAR1Q(nr,nc,mask)
Qphi <- Qs
eigenQphi <- eigen(Qphi, symm=T)
Pphi <- eigenQphi$vector
Dphi <- eigenQphi$value

ns <- nrow(Qphi)
n <- nt * ns
n0 <- n
n1 <- sum(abs(Dtheta)>1e-6)
n2 <- sum(abs(Dxi)>1e-6)
n3 <- sum(abs(Dphi)>1e-6)
n4 <- n2 * n3
N <- n + p + n1 + n2 + n3 + n4
np <- p + 2 * nt + ns + n

getmean <- function(x)
{
  ## mean function
  mu <- as.vector(X %*% x[1:p])
  theta <- x[p + 1:nt]
  xi <- x[ p+ nt + 1:nt]
	phi <- x[ p+ 2*nt + 1:ns]
	delta <- x[p + 2*nt + ns + 1:n]
  mu + rep(theta+xi,each=ns) + rep(phi,nt) + delta
}
loglik <- function(tausq,x)
{
  mm <- getmean(x)
  if(length(id)>0){
    ll <- sum(dnorm(y[-id],mm[-id],1/sqrt(tausq[1]),log=T))
  }
  else{
    ll <- sum(dnorm(y,mm,1/sqrt(tausq[1]),log=T))
  }
  ll
}

genY <- function(tausq,x)
{  mm <- getmean(x)
   sigma <- 1 / sqrt(tausq[1])
   rnorm(length(mm),mm, sigma)
}
genbeta <- function(tausq,x,y)
{
  theta <- x[p + 1:nt]
  xi <- x[ p+ nt + 1:nt]
	phi <- x[ p+ 2*nt + 1:ns]
	delta <- x[p + 2*nt + ns + 1:n]
	r <- y - rep(theta,each=ns) - rep(xi,each=ns) - rep(phi,nt) - delta 
	b <- tausq[1]*as.vector(crossprod(X,r))
	Q <- tausq[1]* crossprod(X)
	genmvnorm2(1,b,Q)
}

genxitheta <- function(tausq,x,y)
{
	beta <- x[1:p]
	phi <- x[p + 2*nt + 1:ns]
	delta <- x[p +2*nt + ns + 1:n]
	mu <- as.vector(X %*% beta)
	r <- y - mu - rep(phi,nt) - delta
	dim(r) <- c(ns , nt)
	b <- apply(r , 2 , sum )
	stopifnot(length(b) == nt)
	Q <- matrix( 0 , 2 * nt , 2 * nt)
	Q[1:nt , 1:nt] <- tausq[1] * ns * diag(nt) + tausq[2] * Qtheta
	Q[1:nt , nt + 1:nt] <- Q[nt + 1:nt , 1:nt] <- tausq[1] * ns * diag(nt)
	Q[nt + 1:nt , nt + 1:nt] <- tausq[1] * ns * diag(nt) + tausq[3] * Qxi
	b <- tausq[1]* c(b , b)
	A <- rep(c(0,1),each=nt)
	e <- 0
	### genmvnorm2 function in m3utils.R 
	### sample multivariate normals with Matrix package
	x <- genmvnorm2(1,b,Q,A,e) 
	c(x)
}

genphi <- function(tausq,x,y)
{
	beta <- x[1:p]; theta <- x[p + 1:nt]
	xi <- x[p + nt + 1:nt]
	delta <- x[p + 2*nt + ns + 1:n]
	mu <- as.vector(X %*% beta)
	dim(mu) <- dim(delta) <- dim(y) <- c(ns,nt)
	b <- apply(y - mu-delta, 1, sum) - (sum(theta) + sum(xi))
	v <- 1 / (tausq[1]*nt + tausq[4] * Dphi)
	m <- crossprod(Pphi,b) * tausq[1] * v
	o <- rnorm(length(m), m , sqrt(v) )
	r <- Pphi %*% o
	r - mean(r)
}

gendelta <- function(tausq, x,y)
{
	beta <- x[1:p]
	theta <- x[p + 1:nt]
	xi <- x[p + nt + 1:nt]
	phi <- x[p + 2*nt + 1:ns]
	mu <- as.vector(X %*% beta)
	eta <- mu + rep(theta + xi , each = ns) + rep(phi , nt)
	b <- y - eta
	v <- 1 / ( tausq[1] + tausq[5]* kronecker(Dxi,Dphi) )
	m <- kmult(t(Pxi),t(Pphi),b) * tausq[1] * v
	o <- rnorm(length(m), m , sqrt(v) )
	r <- kmult(Pxi,Pphi,o)
	dim(r) <- c(ns,nt)
	## print(sum(abs(v)>1e-5))
	c(center(r))	
}

gentausq <- function(x,y)
{
	beta <- x[1:p] 
	theta <- x[p + 1:nt]
	xi <- x[p + nt + 1:nt]
	phi <- x[p + 2*nt + 1:ns]
	delta <- x[p + 2*nt + ns + 1:n]
	mu <- as.vector(X %*%beta)
	R <- y - mu - rep(theta,each=ns) - rep(xi,each=ns) - rep(phi,nt) - delta
	b <- crossprod(R)
	b <- c(b , t(theta) %*% Qtheta %*% theta)
	b <- c(b , t(xi) %*% Qxi %*% xi)
	b <- c(b , t(phi) %*% Qphi %*% phi)
	temp <- kmult(Qxi,Qphi,delta)
	b <- c(b, crossprod(delta,temp))
	b <- parmsb + b / 2
	a <- parmsa + c(n0,n1,n2,n3,n4) / 2
	print(cbind(a,b))
	rgamma(length(a), a , b)
}


ToStore <- rep(0,iters)
Loc <- seq(burnin+1,iters,by=thin)
Ntostore <- length(Loc)
ToStore[Loc] <- 1:Ntostore

Tausq <- array(0,c( Ntostore , 5 , nchains))
dimnames(Tausq) <- list(NULL,paste('tausq',0:4,sep=''),paste('chain',1:nchains))
B <- array(0, c(np , 2 , nchains))
txt <- c(paste('beta',1:p-1,sep=''),paste('theta',1:nt,sep=''),paste('xi',1:nt,sep=''),paste('phi',1:ns,sep=''),paste('delta',1:n,sep=''))
dimnames(B) <- list(txt,c('mean','SSM'),paste('chain',1:nchains))
bid <- c(1:p,sample(p+1:nt,1),sample(p+nt+1:nt,1),
         sample(p+2*nt+1:ns,1),sample(p+2*nt+ns+1:n,2))
Bsub <- array(0, c(Ntostore,length(bid),nchains))
dimnames(Bsub) <- list(NULL,txt[bid],paste('chain',1:nchains,sep=''))
Loglik <- matrix(0, Ntostore, nchains)
DATA <- array(0, c(n, 2,nchains))

if(storeB){
  Bstore <- array(0,c( Ntostore , np , nchains))
}
else{
  Bstore <- NULL
}

if(storeY){
  Ypred <- array(0,c( Ntostore, n, nchains))
}
else{
  Ypred <- NULL
}
## initial values
currB <- matrix(0,nchains,np)
currTausq <- t(inits)
currY <- matrix(0,nchains,n)
for( ch in 1:nchains) {
   currY[ch,] <- y
   currY[ch,id] <- initsy
}

#browser()

## print('start to iterate.')
## ptm0 <- proc.time()
for(j in 1: (iters) ){
	print(j)
	for( k in 1: nchains) {
		## conditioning on missing data
		currB[k,1:p] <- genbeta(currTausq[k,], currB[k,],currY[k,])
		currB[k, p + 1 : (2*nt)] <- genxitheta(currTausq[k,], currB[k,],currY[k,])
		#### plot(ts(cbind(simresult$thetahat, currB[k,1+1:nt])),plot.type='single',ylab='theta',main=j,col=2:1)
		#### plot(ts(cbind(simresult$xi, currB[k,1+nt+1:nt])),plot.type='single',ylab='xi',main=j,col=2:1)
		currB[k,p+2*nt+1:ns] <- genphi(currTausq[k,], currB[k,],currY[k,])
		####plot(ts(cbind(simresult$phi, 	currB[k,1+2*nt+1:ns])),plot.type='single',ylab='phi',main=j,col=2:1)
		currB[k,p+ 2*nt + ns + 1:n] <- gendelta(currTausq[k,], currB[k,], currY[k,])
		#### r <- currB[k,1+2*nt+ns+1:n]
		#### dim(r) <- c(ns,nt)
		####plot(ts(cbind(simresult$delta[1,1,],r[1,])),plot.type='single',ylab='delta11',main=j,col=2:1)
		currTausq[k,] <- gentausq(currB[k,],currY[k,])

		### conditioning on parameters, generate missing data
		tempdata <- genY(currTausq[k,], currB[k,])
		currY[k,id] <- tempdata[id]

		### store
		if( ToStore[j] > 0 ) {
			storeid <- ToStore[j]
			Loglik[storeid,k] <- -2 * loglik(currTausq[k,], currB[k,])
			Tausq[storeid,,k] <- currTausq[k,]
			Bsub[storeid,,k] <- currB[k,bid]
			B[,,k] <- update(B[,,k], currB[k,] , storeid, 0)
			DATA[,,k] <- update(DATA[,,k], tempdata, storeid, 0)
			if(storeB){
			  Bstore[storeid,,k] <- currB[k,]
			}
			if(storeY){
			  Ypred[storeid,,k] <- tempdata
			}
		}
	}
	#ptm <- proc.time() - ptm
	#left <- ptm[3] * (iters + 1 - j )
	#cat(j-1, 'iterations done.' , left,' seconds to go.\n') 
}
Dbar <- apply(Loglik,2,mean)
Tausqbar <- apply(Tausq,c(2,3),mean)
Dhat <- rep(0, nchains)
for( ch in 1:nchains) Dhat[ch] <- -2 * loglik(Tausqbar[,ch],B[,1,ch])
DICout <- list( Dbar = Dbar, Dhat = Dhat, pD = Dbar - Dhat, DIC = 2*Dbar-Dhat)
for( ch in 1:nchains) {
	B[,2,ch] <- B[,2,ch] / (Ntostore - 1)
	DATA[,2,ch] <- DATA[,2,ch] / (Ntostore - 1)
}
timestamp2 <- proc.time()
list(Tausq = Tausq , B = B , Bsub = Bsub,DIC = DICout, Y = DATA, Ypred=Ypred, Ball=Bstore,elapsed =  timestamp2 - timestamp1)

}
