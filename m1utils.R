seasonalQ <- function(per,nt)
{
  
  Q <- matrix(0, nrow=nt, ncol=nt)

  Ks <- matrix( 0, nrow=per, ncol=nt)
  for(i in 1:per ) {
    Ks[i,] <- c( rep(0, i-1), rep(1,per), rep(0, nt-per-i+1) )
  }

  Q[1,]=Ks[1,]
  Q[nt,]=Ks[1,nt:1]
  for(i in 2:per) {
    Q[i,] = apply( Ks[1:i,] , 2 , sum )
    Q[nt-i+1,] = Q[i,nt:1]
  }
  for(i in (per+1):(nt-per)){
    Q[i,] = c(rep(0,i-per),c(1:per,(per-1):1),rep(0,nt-per+1-i))
  }
  Q
}

SRW1Q <- function(per,nt)
{
	K <- matrix(0, nrow=nt-per,ncol=nt)
	for(i in 1:(nt - per) ) K[i, c(i,i+per)] <- c( -1,1)
	Q <- crossprod(K)
	Q
}


RW1Q <- function(n)
{
	K <- matrix(0, nrow = n -1 , ncol=n)
	for( i in 1:nrow(K))
		K[i,c(i,i+1)] <- c(1,-1)
	Q <- crossprod(K)
	Q
}


RW2Q <- function(n)
{
  Q <- matrix(0, nrow=n, ncol=n)
  Q[1,1:3] <- Q[n, n:(n-2)] <- c(1,-2,1)
  Q[2,1:4] <- Q[n-1, n:(n-3)] <- c(-2,5,-4,1)
  for(j in 3:(n-2)) 
    Q[j,] <- c(rep(0,j-3), c(1,-4,6,-4,1), rep(0, n-j-2))
  
  print("Q")
  print(Q)
  
}

CAR1Q <- function(nr,nc,mask=NULL,nb=NULL,type='rook')
{
	## construct spatial component
	##  nr*nc x nr*nc precision matrix
	##  mask from zero
	require(spdep)
	if(is.null(nb)){
		nb1111 <- cell2nb(nr,nc,type=type)
	}else{
		nb1111 <- nb
	}
	if(!is.null(mask)){
		todrop <- as.logical(c(mask))
		nb1111 <- droplinks(nb1111,todrop)
	}
	mat1111 <- nb2mat(nb1111,style="B",zero.policy=T)
	temp <- apply(mat1111,2,sum)
	diag(mat1111) <- -temp
	Qs <- -mat1111
	if(!is.null(mask)) Qs <- Qs[!todrop,!todrop]
	rm(nb1111)
	rm(mat1111)
	rm(temp)
	Qs
}

newseasonalQ <- function(per)
{
	### Cowles new seasonal model
  half <- per %/% 2
  K <- matrix( 0, nrow = per, ncol=per)
  for( p in 1:( half -1) ) {
       K[p, p+1] <- 1
       K[p, per-p+1] <- -1
  }
  
  if(!( per %% 2 )) { # even number
      K[ half, 1] <- K[half, half+1] <- 1
  }
  else {
        K[ half, half+1] <- 1
        K[half, half+2] <- -1
  }

  Q <- t(K) %*% K
  Q
}

KRW1Q <- function(n)
{
	K <- matrix(0, nrow = n -1 , ncol=n)
	for( i in 1:nrow(K))
		K[i,c(i,i+1)] <- c(1,-1)
	K
}

KseasonalQ <- function(per,n)
{
	# for seasonal trend
  Ks <- matrix( 0, nrow=n-per+1, ncol=n)
  for(i in 1:(n-per+1) ) Ks[i,] <- c( rep(0, i-1), rep(1,per), rep(0, n-per-i+1) )
	Ks
}

KnewseasonalQ <- function(per)
{
	### Cowles new seasonal model
  half <- per %/% 2
  K <- matrix( 0, nrow = half, ncol=per)
  for( p in 1:( half -1) ) {
       K[p, p+1] <- 1
       K[p, per-p+1] <- -1
  }
  
  if(!( per %% 2 )) { # even number
      K[ half, 1] <- K[half, half+1] <- 1
  }
  else {
        K[ half, half+1] <- 1
        K[half, half+2] <- -1
  }
	K
}

KSRW1Q <- function(per,nt)
{
	K <- matrix(0, nrow=nt-per,ncol=nt)
	for(i in 1:(nt - per) ) K[i, c(i,i+per)] <- c( -1,1)
	K
}

KRW2Q <- function(n)
{
# construct design matrices
# for RW2
  Kt <- matrix( 0, nrow=n-2, ncol=n)
  for( i in 1:(n-2) )  Kt[i,] <- c(rep(0,i-1), c(1,-2,1), rep(0, n - 2 - i))
	Kt
}

KCAR1Q <- function(nr,nc,mask=NULL,type='rook')
{
	Qs <- CAR1Q(nr,nc,mask=mask,type=type)
	eigenQs <- eigen(Qs,symm=T)
	lambda <- round(eigenQs$value,10)
	idx <- which(lambda >0)
	lambda <- lambda[idx]
	Ps <- eigenQs$vector[,idx]
	Ks <- Ps %*% diag(sqrt(lambda))
	##plot(Ks %*% t(Ks) - Qs)
	t(Ks)
}

makeSseasonalRue <- function(per, nt)
{
	m <- matrix(0, nrow=nt, ncol=( per - 1) )
	for( i in 1 : nt )
	{
		j <- i %% per
		if( j == 0 ) m[i,] <- -1
		else m[i , j ] <- 1
	}
	m
}

makeSseasonalRW <- function(per,nt)
{
	m <- matrix(0, nrow=nt, ncol = per)
	for( i in 1:nt )
	{
		j <- (i - 1) %% per + 1
		m[i, j ] <- 1
	}
	m
}

SnewseasonalQ <- function(per)
{
	### Cowles new seasonal model
  half <- per %/% 2
  K <- matrix( 0, nrow = half, ncol=per)
  for( p in 1:( half -1) ) {
       K[p, p+1] <- 1
       K[p, per-p+1] <- 1
  }
  
  if(!( per %% 2 )) { # even number
      K[ half, 1] <- 1
      K[half, half+1] <- -1
  }
  else {
        K[ half, half+1] <- 1
        K[half, half+2] <- 1
  }
	t(K)
}


simu <- function(k,Q, title){
  s = rep(0,nrow(Q))
  spd = eigen(k*Q)
  L = spd$values
  V = spd$vectors
  for(j in 1:length(L)){
    if( abs(L[j]) > 1e-10 ) {
      s = s + rnorm(1 , 0 , 1/sqrt(L[j]) ) * V[,j]
    }
  }
  plot(ts(s),main=title)
  s
}

nb2WB2 <- function(nbr,todrop = NULL)
{
if(!is.null(todrop)){
	nbr <- droplinks(nbr,todrop)
	table0 <- which(!todrop)
	adj <- unlist(nbr)
	adj <- adj[adj>0]
	adj <- match(adj, table0, nomatch=0)
	num <- unlist(lapply(nbr, length))
	num <- num[!todrop]
	sumnumneighs <- sum(num)
	weights <- rep(1 , sumnumneighs)
	l <- list(adj = adj, weights = weights, num = num)
}else 
	l <- nb2WB(nbr)
	l
}