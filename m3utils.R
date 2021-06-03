makeQseasonalRue <- function(per,nt)
{
	# called by main func
	K <- matrix( 0, nrow=nt-per+1, ncol=nt )
	for(i in 1:nrow(K) )
        K[i, i:(i+per-1) ] <- rep(1,per)
	Q <- t(K) %*% K
Q
}


makeRW2seasonalRuewierd4 <-  function(per, nt)
{
	# start with Rue-type seasonal Q
	Qseas <-makeQseasonalRue(per,nt)

	# make it invertible
	Qseasposdef <- Qseas
	diag(Qseasposdef) <- diag(Qseasposdef) + 1

	Sigma <- solve(Qseasposdef)

	# compute conditional cov matrix given Bx = 0
	# B has vectors with 0 eigenvals in both seasonal & RW2

	B <- t(cbind( makeS(per,nt), rep(1,nt), 1:nt))

	Sigmasuboff <- t(B) %*% solve( (B) %*% Sigma %*% t(B))%*% B

	Sigmaconst <- Sigma - Sigma %*% Sigmasuboff %*% Sigma


	# now get equivalent Q mat
	eigenSc <- eigen(Sigmaconst)

	evals <- (eigenSc$values)

	k <- length(evals[abs(evals) < 0.0000000001])

	evalsinv <- c( 1/evals[1:(nt-k)], rep(0,k) )

	Q <- eigenSc$vectors %*% diag(evalsinv) %*% t(eigenSc$vectors)

	list(Sigmaconst=Sigmaconst, Q=Q)
}
	
	fakeinvert <- function(mat) {
		# mat has to be symmetric nonneg def
		nt <- nrow(mat)
		eigenmat <- eigen(mat)
		evals <- eigenmat$values
		k <- length( evals[ abs(evals) < 0.0000000001] )
		evalsinv <- c(1/evals[1:(nt-k)],rep(0,k) )
		fakeinverse <- eigenmat$vectors %*% diag(evalsinv) %*%t(eigenmat$vectors)
		fakeinverse
	}

makeRW2seasonalRuewierd5 <- function(per, nt)
{
	Qseas <-makeQseasonalRue(per,nt)
	Sigma <- fakeinvert(Qseas)
	B <- t(cbind( makeSseasonalRue(per,nt), rep(1,nt), 1:nt))
	temp <- B %*% Sigma %*% t(B)
	tempinv <- fakeinvert(temp)
	Sigmasuboff <- t(B) %*% tempinv %*% B
	Sigmaconst <- Sigma - Sigma %*% Sigmasuboff %*% Sigma
	Q <- fakeinvert(Sigmaconst)
	
	list(Sigmaconst=Sigmaconst, Q=Q)
}

runif.simplex <- function(n) {
   ## draws uniformly from standard simplex with nvert vertices
   ## require(MCMCpack)
   rdirichlet(1, rep(1, n))
}

rdirichlet <- function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

makeFirstSimplex <- function(kappa, mult)
## kappa:   current kappa vector
## mult :   the multiplier to get the size of the random simplex
{
   vertices <- diag(length(kappa))

   if (mult < 1) {
      ## get the coordinates of the simplex shrinked towards vertex 1
      vertices[,-1] <- vertices[,-1] + (1 - mult) * (vertices[,1] - vertices[,-1])

      ## first, generate a point uniformly in appropriate size triangle
      ## then translate vertices to make x correspond to random point
      rbaryc <- runif.simplex(length(kappa))
      vertices <- vertices + kappa - as.vector(vertices %*% t(rbaryc))
   }

   vertices
}

shrinkSimplex <- function(bx, bc, cx, cc, vertices)
## bx: coefficients of vertices to produce x (barycentric coordinates of x)
## bc: coefficients of vertices to produce cand (barycentric coords of cand)
## cx, cc:  Cartesian coords of x and cand
## vertices: vertices of the current simplex
{
   for (i in seq(bc)[bc < bx]) {
      vertices[,-i] <- vertices[,-i] + bc[i] * (vertices[,i] - vertices[,-i])
      bc <- solve(vertices, cc)
   }

   vertices
}

sliceSimplex <- function(theta, mpdfun, log=TRUE, f.theta, mult, ...)
## theta:    current point
## mpdfun:   marginalized posterior density up to a scale or additive constant
## log:      work on log fun or not
## f.theta:  mpdfun(theta) or NA to avoid using that from previous iteration
{
   if(is.na(f.theta)) {evalresult<- mpdfun(kappa=theta,...); f.theta <-evalresult$value}
   y <- ifelse(log, f.theta - rexp(1), runif(1, 0, f.theta))
	 #print(c(f.theta,y))	 
   ## simplex for generating kappa
   kappa.old <- theta
   vertices <- makeFirstSimplex(kappa.old, mult)
   kappa.old.b <- solve(vertices, kappa.old)
   newevals <- 0
   repeat {
      ## sample kappa.cand
      kappa.cand.b <- runif.simplex(length(kappa.old))
      kappa.cand <- as.vector(vertices %*% t(kappa.cand.b))

      # out of boundary indicator
      ob.kappa <- min(kappa.cand) < 0 || max(kappa.cand) > 1

      if (!ob.kappa) { ## cand is in the boundary
         cand <- kappa.cand
         evalresult <- mpdfun(kappa=cand,...)
         if (evalresult$value > -Inf) {  ## valid x0 to mpdfun
            newevals <- newevals + 1
         }
         # print(rbind(c(theta,y),c(cand,evalresult$value)))
				 ## here is the only exit of the loop
         if (y < evalresult$value) {
            evalresult$theta <- cand
            evalresult$newevals <- newevals
            return(evalresult)
         }
      }

      vertices <- shrinkSimplex(kappa.old.b, kappa.cand.b, kappa.old,
                                kappa.cand, vertices)
      kappa.old.b <- solve(vertices, kappa.old)
   }
}

kmult.0 <- function( A , B , z )
{
	#### (A \otimes B) %*% z
	#### A  --- (nt x nt ) time component
	#### B  --- (ns x ns ) spatial component
	#### z  --- (ns * nt ) x 1 vector , ordered site within time
	#### version 0
	nt <- nrow(A)
	ns <- nrow(B)
	stopifnot(nrow(A)==ncol(A))
	stopifnot(nrow(B)==ncol(B))
	Z <- matrix(c(z), nrow = nt, byrow = T )
	Bt <- t(B)
	y <- rep( 0, length(z) )
	for( i in 1 : nt ) {
		temp <- crossprod(Z, c(A[i,]))
		y[ (i-1)*ns + 1:ns ] <- as.vector(crossprod(Bt , temp ))
	}
	y
}
kmult <- function( A , B , x )
{
	#### (A \otimes B) %*% x
	#### A  --- (ntheta x nt ) time component
	#### B  --- (nphi x ns ) spatial component
	#### z  --- (ns * nt ) x 1 vector , ordered site within time
	#### version 1 
	#### due to Jain, Anil K. (1989), Fundamentals of Digital Image Processing
	####        (A \otimes B)%*% x = B X A^T  
	#### where x is row ordered into column vectors
	X <- matrix(c(x), nrow = ncol(B))
	Z <- B %*% X
	Y <- tcrossprod(Z,A)
	c(as.matrix(Y)) 	
}

kmult.back <- function( A , B , x )
{
	#### (A \otimes B) %*% x
	#### A  --- (ntheta x nt ) time component
	#### B  --- (nphi x ns ) spatial component
	#### z  --- (ns * nt ) x 1 vector , ordered site within time
	#### version 1 
	#### due to Jain, Anil K. (1989), Fundamentals of Digital Image Processing
	####        (A \otimes B)%*% x = B X A^T  
	#### where x is row ordered into column vectors
	require(Matrix)
	A <- Matrix(A)
	B <- Matrix(B)
	X <- matrix(c(x), nrow = ncol(B))
	Z <- B %*% X
	Y <- tcrossprod(Z,A)
	c(as.matrix(Y)) 	
}

update <- function(muSS, x,iter,burnin)
{
	## Algorithm due to Knuth
	## update sample mean and sum of square deviation from mean in sequence
	stopifnot(!is.null(dim(muSS)))
	stopifnot(nrow(muSS) == length(x))
	#browser()
	mu <- muSS[,1]
	SSd <- muSS[,2]
	m <- iter - burnin
	if( m > 0 )
	{	
		delta <- x - mu
		mu <- mu + delta / m
  	SSd <- SSd + delta * (x - mu)
	}
	cbind(mu, SSd)
}

simuinter <- function(Qt , Qs , tau)
{
  ## simulate interaction fields based on Qt \otimes Qs structure matrix
  ## Algorithm due to Rue Algorithm 3.1
  ## Qt  : temporal structure matrix
  ## Qs  : spatil structure matrix
  ## tau : true precision
  
  nt <- nrow(Qt); ns <- nrow(Qs)
  n <- nt * ns
  eigenQt <- eigen(Qt, symm = T)
  Pt <- eigenQt$vector; Dt <- eigenQt$value
  eigenQs <- eigen(Qs, symm = T)
  Ps <- eigenQs$vector; Ds <- eigenQs$value
  
  Delta <- round(kronecker(Dt, Ds),10)    
  Delta <- tau * Delta
  ## temp <- matrix(0,nrow=n,ncol=n)
  delta <- rep(0, n) 
	for( j in 1:n ){
		ct <- ceiling ( j / ns )
		cs <- (j -1) %% ns + 1
		print(c(ct,cs))
		stopifnot( (cs + (ct-1)*ns) == j )
		v <- kronecker(Pt[,ct], Ps[,cs])
		if( Delta[j] > 0) {
			#temp <- temp + Delta[j] * tcrossprod(v)
			delta = delta + rnorm(1 , 0 , 1/sqrt(Delta[j]) ) * v
		}
	}
  delta
}	
  	
trend <- function(Qt,Qs, y)
{
	## Qt
	## Qs
	## y
	## consistent trend 
	nt <- nrow(Qt); ns <- nrow(Qs); n <- nt *ns
	eigenQt <- eigen(Qt, symm = T)
  Pt <- eigenQt$vector; Dt <- eigenQt$value
  eigenQs <- eigen(Qs, symm = T)
  Ps <- eigenQs$vector; Ds <- eigenQs$value
	Delta <- round(kronecker(Dt, Ds),10)
	id <- which(Delta==0)
	
	Ny <- rep(0,n)
	for( j in id ){
		ca <- ceiling ( j / ns )
		cb <- (j -1) %% ns + 1
		print(c(ca,cb))
		stopifnot( (cb + (ca-1)*ns) == j )
		v <- kronecker(Pt[,ca],Ps[,cb])
		Ny = Ny + c(crossprod(v,y)) * v
	}
	Ny
 ## i didn't change 1x1 matrix to numeric	
 ## i used non-zero eigenvalues
}

center <- function(m)
{
	mm <- scale(m,scale=F)
	mm <- t(mm)
	mm <- t(scale(mm,scale=F))
	mm
}	

genmvnorm1 <- function(n,b,Q)
{
	### naive sample from multivariate normal distribution
	### n     :  number of samples
	### b , Q : defined by p(x) \propto \exp( -.5 x^T Q x + b^T x)
	require(mvtnorm)
	require(Matrix)
	Q <- Matrix(Q)
	Sigma <- solve(Q)
	mu <- Sigma %*% b
	x <- rmvnorm(n, as.vector(mu), as.matrix(Sigma) )
	x
}
genmvnorm2new <- function(n,b,Q,A=NULL,e=NULL)
{
	### Fast sample from multivariate normal
	### Algorithm due to Rue 2001 JRSSB
	###   n   : number of samples
	### b , Q : defined by p(x) \propto \exp( -.5 x^T Q x + b^T x)
	### Q     : usually sparse
	### A,  e : defining constraints Ax = e
      ### simplifed based on spam Aug 20 2009 
	ch <- chol(Q)
      mu <- backsolve(ch,forwardsolve(ch,b))
	eps <- matrix( rnorm(n* ncol(Q)) , ncol =n )
	x <- t(backsolve(ch , eps ))
	x <- sweep(x , 2 , mu , '+')
      if( !is.null(A) && !is.null(e) )
      {
  	  if(is.null(dim(A))) dim(A) <- c(1,ncol(Q))
  	  V <- backsolve(ch, forwardsolve(ch,t(A)))
        W <- as.matrix(A %*% V)
  	  U <- solve(W , t(V) )
  	  C <- t(apply(x , 1 , function(v) 
              as.vector(as.matrix(crossprod( U, A %*% v - e ) ))
              ))
  	   x <- x - C
      }
      x
}

genmvnorm2 <- function(n,b,Q,A=NULL,e=NULL)
{
	### Fast sample from multivariate normal
	### Algorithm due to Rue 2001 JRSSB
	###   n   : number of samples
	### b , Q : defined by p(x) \propto \exp( -.5 x^T Q x + b^T x)
	### Q     : usually sparse
	### A,  e : defining constraints Ax = e
	require(Matrix)
	Q <- Matrix(Q) 
	ch <- chol(Q, pivot=T)
	oo <- attr(ch, "pivot")
	if( !is.null(oo) ) {
		o <- order(oo)
		mu <- backsolve(ch , forwardsolve(ch, b[oo] , upper=T, trans=T) )
	} else
	{
		mu <- backsolve(ch , forwardsolve(ch, b , upper=T, trans=T) )
	}		
	eps <- matrix( rnorm(n* ncol(Q)) , ncol =n )
	x <- t(backsolve(ch , eps ))
  x <- sweep(x , 2 , mu , '+')
  dx <- dim(x)
  if( !is.null(oo) ) x <- x[,o]
  dim(x) <- dx
  if( !is.null(A) && !is.null(e) )
  {
  	if(is.null(dim(A))) dim(A) <- c(1,ncol(Q))
  	if( !is.null(oo) ){
  		Ap <- t(A)[oo,]
  		if(is.null(dim(Ap))) dim(Ap) <- c( ncol(Q) , 1 )
  		Vp <- backsolve(ch, forwardsolve(ch, Ap, upper=T, trans=T) )
  		V <- Vp[o,]
		} else
		{
			V <- backsolve(ch, forwardsolve(ch, t(A) , upper = T , trans =T) )
		}	
  	      W <- A %*% V
  		U <- solve(W , t(V) )
  	C <- t(apply(x , 1 , function(v) 
              as.vector(as.matrix(crossprod( U, A %*% v - e ) ))
           ))
  	x <- x - C
  }	
  x
}

genmvnorm3 <- function(n,b,Q,oo,A=NULL,e=NULL)
{
	### Fast sample from multivariate normal resuing permutation
	### Algorithm due to Rue 2001 JRSSB
	###   n   : number of samples
	### b , Q : defined by p(x) \propto \exp( -.5 x^T Q x + b^T x)
	### Q     : usually sparse
	### oo    : permutation of nodes to reduce fill-in
	### A,  e : defining constraints Ax = e
	require(Matrix)
	Q <- Matrix(Q[oo,oo])
	ch <- chol(Q)
	o <- order(oo)
	mu <- backsolve(ch , forwardsolve(ch, b[oo] , upper=T, trans=T) )
	eps <- matrix( rnorm(n* ncol(Q)) , ncol =n )
	x <- t(backsolve(ch , eps ))
  x <- sweep(x , 2 , mu , '+')
  dx <- dim(x)
  x <- x[,o]
  dim(x) <- dx
  if( !is.null(A) && !is.null(e) )
  {
  	if(is.null(dim(A))) dim(A) <- c(1,ncol(Q))
  	Ap <- t(A)[oo,]
  	if(is.null(dim(Ap))) dim(Ap) <- c( ncol(Q) , 1 )
  	Vp <- backsolve(ch, forwardsolve(ch, Ap, upper=T, trans=T) )
  	V <- Vp[o,]
  	W <- A %*% V
  	U <- solve(W , t(V) )
  	C <- t(apply(x , 1 , function(v) crossprod( U, A %*% v - e ) ))
  	x <- x - C
  }	
  x
}

genmvnorm3b <- function(n,Y, X,gamma,oo,A=NULL,e=NULL)
{
	### Fast sample from multivariate normal resuing permutation
	### Algorithm due to Rue 2001 JRSSB
	###   n   : number of samples
	### X ,Y  : defined by Hodge reparameterization Y = X \Theta + E
	### gamma : diagonal entris of Gamma defined by Hodge repara. E ~ N(0,Gamma) 
	### oo    : permutation of nodes to reduce fill-in
	### A,  e : defining constraints Ax = e
	require(Matrix)
	b <- crossprod(X, gamma * Y )
	k <- Diagonal(x = gamma)
	Q <- t(X) %*% k %*% X  + 1e-7 * Diagonal(ncol(X))
	Q <- Matrix(Q[oo,oo])
	ch <- chol(Q)
	o <- order(oo)
	mu <- backsolve(ch , forwardsolve(ch, b[oo] , upper=T, trans=T) )
	eps <- matrix( rnorm(n* ncol(Q)) , ncol =n )
	x <- t(backsolve(ch , eps ))
  x <- sweep(x , 2 , mu , '+')
  dx <- dim(x)
  x <- x[,o]
  dim(x) <- dx
  if( !is.null(A) && !is.null(e) )
  {
  	if(is.null(dim(A))) dim(A) <- c(1,ncol(Q))
  	Ap <- t(A)[oo,]
  	if(is.null(dim(Ap))) dim(Ap) <- c( ncol(Q) , 1 )
  	Vp <- backsolve(ch, forwardsolve(ch, Ap, upper=T, trans=T) )
  	V <- Vp[o,]
  	W <- A %*% V
  	U <- solve(W , t(V) )
  	C <- t(apply(x , 1 , function(v) crossprod( U, A %*% v - e ) ))
  	x <- x - C
  }	
  x
}

cg <- function(A,b,x,imax,tol)
{
  ## conjudate gradient method 
  ## algorithm due to Jonathan Richard Shewchuk
	## A      : 
	## b      :  define A x  = b
	## x      :  starting value
	## imax   :  maximum number of iterations
	## tol    :  error tolerance
	## Output :  solution to A x = b
	i <- 0
	r <- b - A %*% x
	d <- r
	delta.new <- c(as.matrix(crossprod(r)))
	delta0 <- delta.new
	while( (i < imax) && ( delta.new > tol^2 ) )
	{
		q <- A %*% d
		alpha <- delta.new / c(as.matrix(crossprod(d , q)))
		x <- x + alpha * d
		if( (i %% 50) == 0)
		{
			r <- b - A %*% x
		}else
		{
			r <- r - alpha * q
		}
		delta.old <- delta.new
		delta.new <- c(as.matrix(crossprod(r)))
		print(delta.new)
		beta <- delta.new / delta.old
		d <- r + beta * d
		i <- i + 1
	}
	c(as.matrix(x))
}
	
genmvnorm4 <- function(X, Gamma, Y, x , imax, tol)
{
	### High dimensional multivariate normal sampling
	### algorithm due to Wikle et al 01 JASA 
	### Input
	### X      : N x n SMCMC design matrix
	### Gamma  : vector of length N, diagonal of precisions
	### Y      : vector of length N, data
	### x      : inital value for Conjugate gradient method
	### imax   : maximum number of iterations
	### tol    : numerical tolerance for CG method
	### output
	###    a sample from N( (X^t Gamma X)^{-1} X^t Gamma Y, X^t Gamma X)
	e <- rnorm(length(Gamma))
	f <- crossprod(X, sqrt(Gamma) * e)
	g <- crossprod(X, Gamma * Y )
	k <- Diagonal(x = Gamma)
	Q <- t(X) %*% k %*% X
	out <- cg(A=Q,b=g+f,x=x,imax=imax,tol=tol)
}

geninter <- function(b,Qa,Qb,Ka,Kb,ky,k,D,x,imax,epsilon=0.0005)
{ 
  ### simulate from 
  ### p(x) \propto \exp \left( -\frac{1}{2} xtQx + btx\right)
  ### with Q \equiv ky * D + k * Qa \otimes Qb
  ### b       : defined above(don't miss constants)
  ### Ka      : defined via Qa \equiv crossprod(Ka)
  ### Kb      : devind via Qb \equiv crossprod(Kb)
  ### ky      : measurement error precisions
  ### k       : random effect precisions
  ### D       : Diagonal entries of diagonal matrix D
	### x       : inital value for Conjugate gradient method
	### imax    : maximum number of iterations
	### epsilon : numerical tolerance for CG method
	kmult2 <- function (k,A,B,ky,D,x)
	{ 
		#### k * A \otimes B x + ky * D * x
		k * kmult(A,B,x) + ky * D * x
	}
	e1 <- rnorm(length(D))
	e2 <- rnorm(nrow(Ka) * nrow(Kb))
	f <- sqrt(ky*D)*e1 + kmult(t(Ka),t(Kb),sqrt(k)*e2)
  b <- b + f
  tol <- epsilon * crossprod(b)
  ####################################################
  ### use Conjugate gradient algorithm to solve
  ### Q x = b
  i <- 0
	#r <- b - A %*% x
	r <- b - kmult2(k,Qa,Qb,ky,D,x)
	d <- r
	delta.new <- c(as.matrix(crossprod(r)))
	delta0 <- delta.new
	while( (i < imax) && ( delta.new > tol ) )
	{
		#q <- A %*% d
		q <- kmult2(k,Qa,Qb,ky,D,d)
		alpha <- delta.new / c(as.matrix(crossprod(d , q)))
		x <- x + alpha * d
		if( (i %% 50) == 0)
		{
			#r <- b - A %*% x
			r <- b - kmult2(k,Qa,Qb,ky,D,x)
		}else
		{
			r <- r - alpha * q
		}
		delta.old <- delta.new
		delta.new <- c(as.matrix(crossprod(r)))
		## print(delta.new)
		beta <- delta.new / delta.old
		d <- r + beta * d
		i <- i + 1
	}
	c(as.matrix(x))
}

geninter2 <- function(b,Qa,Qb,Ka,Kb,k,D,x,imax,epsilon=0.0005)
{ 
  ### simulate from 
  ### p(x) \propto \exp \left( -\frac{1}{2} xtQx + btx\right)
  ### with Q \equiv D + k * Qa \otimes Qb
  ### b       : defined above(don't miss constants)
  ### Ka      : defined via Qa \equiv crossprod(Ka)
  ### Kb      : devind via Qb \equiv crossprod(Kb)
  ### k       : random effects precisions
  ### D       : Diagonal entries of diagonal matrix D
	### x       : inital value for Conjugate gradient method
	### imax    : maximum number of iterations
	### epsilon : numerical tolerance for CG method
	kmult2 <- function (k,A,B,D,x)
	{ 
		#### A \otimes B x + D * x
		k * kmult(A,B,x) + D * x
	}
	e1 <- rnorm(length(D))
	e2 <- rnorm(nrow(Ka) * nrow(Kb))
	f <- sqrt(D)*e1 + kmult(t(Ka),t(Kb),sqrt(k) *e2)
  b <- b + f
  tol <- epsilon * crossprod(b)
  ####################################################
  ### use Conjugate gradient algorithm to solve
  ### Q x = b
  i <- 0
	r <- b - kmult2(k,Qa,Qb,D,x)
	d <- r
	delta.new <- c(as.matrix(crossprod(r)))
	delta0 <- delta.new
	while( (i < imax) && ( delta.new > tol ) )
	{
		q <- kmult2(k,Qa,Qb,D,d)
		alpha <- delta.new / c(as.matrix(crossprod(d , q)))
		x <- x + alpha * d
		if( (i %% 50) == 0)
		{
			r <- b - kmult2(k,Qa,Qb,D,x)
		}else
		{
			r <- r - alpha * q
		}
		delta.old <- delta.new
		delta.new <- c(as.matrix(crossprod(r)))
		## print(delta.new)
		beta <- delta.new / delta.old
		d <- r + beta * d
		i <- i + 1
	}
	c(as.matrix(x))
}

geninter3 <- function(b,Va,Vb,Lambda,W,tau0,taux,x,imax,epsilon=0.0005)
{ 
  ### simulate from 
  ### p(x) \propto \exp \left( -\frac{1}{2} xTQx + btx\right)
  ### with Q \equiv taux Lambda + tau0 V^T W V
  ###      b \equiv tau0 y^T W V
  ###      V \equiv Va \otimes Vb  
  ### input:
  ### b       : defined above(don't miss constants)
  ### Va      : defined via  
  ### Vb      : devindd via  spectral decomposition of Qa and Qb
  ### Lambda  : defined via
  ### W       : diagonal entries of weight matrix W 
  ### tau0    : measurement error precisions
  ### taux     : random effect precisions
  ### x       : inital value for Conjugate gradient method
  ### imax    : maximum number of iterations
  ### epsilon : numerical tolerance for CG method
	
  Qmult <- function (x,Va,Vb,Lambda,W,tau0,taux)
  { 
    tmp <- W * kmult(Va,Vb,x)
    taux * Lambda * x + tau0 * kmult(t(Va), t(Vb), tmp)
  }
  e1 <- sqrt(taux * Lambda) * rnorm(length(Lambda))
  e2tilde <- rnorm(nrow(Va) * nrow(Vb),sd = sqrt(tau0*W))
  e2 <- kmult(t(Va),t(Vb), e2tilde)
  e <- e1 + e2 
  b <- b + e
  tol <- epsilon * crossprod(b)
  ####################################################
  ### use Conjugate gradient algorithm to solve
  ### Q x = b
  i <- 0
  r <- b - Qmult(x,Va,Vb,Lambda,W,tau0,taux)
  d <- r
  delta.new <- c(as.matrix(crossprod(r)))
  delta0 <- delta.new
  while( (i < imax) && ( delta.new > tol ) )
  {
	q <- Qmult(d,Va,Vb,Lambda,W,tau0,taux)
	alpha <- delta.new / c(as.matrix(crossprod(d , q)))
	x <- x + alpha * d
	if( (i %% 50) == 0)
	{
		r <- b - Qmult(x,Va,Vb,Lambda,W,tau0,taux)
	}else
	{
		r <- r - alpha * q
	}
	delta.old <- delta.new
	delta.new <- c(as.matrix(crossprod(r)))
	beta <- delta.new / delta.old
	d <- r + beta * d
	i <- i + 1
   }
  c(as.matrix(x))
}

geninter4 <- function(b,Va,Vb,Lambda,W,tau0,taux,x,imax,epsilon=0.0005)
{ 
  ### simulate from 
  ### p(x) \propto \exp \left( -\frac{1}{2} xTQx + btx\right)
  ### with Q \equiv tau0 W + taux V^T Lambda V
  ###      Lambda \equiv Lambda_a \otimes Lambda_b 
  ###      V \equiv Va \otimes Vb  
  ### input:
  ### b       : defined above(don't miss constants)
  ### Va      : defined via  
  ### Vb      : devindd via  spectral decomposition of Qa and Qb
  ### Lambda  : defined via
  ### W       : diagonal entries of weight matrix W 
  ### tau0    : measurement error precisions
  ### taux     : random effect precisions
  ### x       : inital value for Conjugate gradient method
  ### imax    : maximum number of iterations
  ### epsilon : numerical tolerance for CG method
  Qmult <- function (x,Va,Vb,Lambda,W,tau0,taux)
  { 
    tmp <- Lambda * kmult(t(Va),t(Vb),x)
    tau0 * W * x + taux * kmult(Va,Vb,tmp) 
  }
  e1 <- sqrt(tau0 * W) * rnorm(length(W))
  e2tilde <- rnorm(nrow(Va) * nrow(Vb),sd = sqrt(taux* Lambda))
  e2 <- kmult(Va,Vb, e2tilde)
  e <- e1 + e2 
  b <- b + e
  tol <- epsilon * crossprod(b)
  ####################################################
  ### use Conjugate gradient algorithm to solve
  ### Q x = b
  i <- 0
  r <- b - Qmult(x,Va,Vb,Lambda,W,tau0,taux)
  d <- r
  delta.new <- c(as.matrix(crossprod(r)))
  delta0 <- delta.new
  while( (i < imax) && ( delta.new > tol ) )
  {
	q <- Qmult(d,Va,Vb,Lambda,W,tau0,taux)
	alpha <- delta.new / c(as.matrix(crossprod(d , q)))
	x <- x + alpha * d
	if( (i %% 50) == 0)
	{
		r <- b - Qmult(x,Va,Vb,Lambda,W,tau0,taux)
	}else
	{
		r <- r - alpha * q
	}
	delta.old <- delta.new
	delta.new <- c(as.matrix(crossprod(r)))
	beta <- delta.new / delta.old
	d <- r + beta * d
	i <- i + 1
   }
  c(as.matrix(x))
}

gengamma <- function(n,alpha,beta=1,tau)
{
  ## generation of truncated gamma random variable
  stopifnot( alpha >=1)
  if(alpha==1) {
  	 Lstar <- tau * beta
       out <- rexp(n) + Lstar
	 res <- list(sample= out / beta , n = n)
  }else{
	res <- gengamma1(n,alpha,beta,tau)
  }
  res
}

gengamma1 <- function(n,alpha,beta=1,tau)
{
  ## generation of truncated gamma random variable
  ## Algorithm due to Dagpunar 79
  ## n   : number of random variates
  ## alpha   : 
  ## beta   : defined by p(x) \propto x^{alpha-1} \exp{-beta x}
  ## tau   :                        I_{tau,\infty}(x)
  ## setup
  i <- 1
  j <- 0	
  mode <- alpha - 1
  Lstar <- tau * beta
  out <- rep(0,n)
  if ( mode > Lstar ) {
     repeat{
        y <- rgamma(1,alpha,1)
	  j <- j + 1
	  if( y > Lstar ){
		out[i] <- y
		i <- i + 1
	  }
	  if( i > n) break
     }	
   } 
   else
   {
 	 lambda <- (Lstar - alpha + sqrt( (Lstar-alpha)^2 + 4 * Lstar )) / (2*Lstar)
       repeat{
	      u <- runif(2)
		y <- - log(u[1]) / lambda + Lstar
            j <- j + 1
            lacct <- log(u[2]) + (1 - lambda) * y + (alpha-1)*( log( alpha -1) - log(1-lambda) - 1 - log(y)) 
		if( lacct <=0 ) {
			out[i] <- y
			i <- i + 1
		}
		if( i > n) break
	   }
    }
 list(sample= out / beta , n = j)
}
gengammaold <- function(n,a,b,A)
{
  ## generation of truncated gamma random variable
  ## Algorithm due to Geweke 91
  ## n   : number of random variates
  ## a   : 
  ## b   : defined by p(x) \propto x^{a-1} \exp{-b x}
  ## A   :                              I_{A^(-2),\infty}(x)
  #### bug a/b big, r overflow to infty
  ###  fix ? calculate logr instead   
	alpha <- (a -1) / b
	out <- rep(0,n)
	if ( A^(-2) < alpha) {
		source = 'gamma'
	} 
	else
	{
		source = 'exp'
		if( A^(-2) * b > 1 ) {
			lambda <- A^2
			if( A^(-2) * b <= a ){
			  temp <- (a-1)/(b-A^{2})
			  #r <- temp^(a-1) * exp(-(b-A^(2))*temp)
			logr <- (a-1)*log(temp) - (b - A^2)*temp
			} else{
			  #r <- A^(-2*(a-1)) * exp(-(A^(-2)*b - 1))
			logr <- -2*a - A^2*b + 3
			}
		}else{
			lambda <- b / 2
			temp <- 2*(a-1)/b
			if( temp >= A^(-2) ){
				#r <- temp^(a-1) * exp(-(a-1))
			logr <- (a-1)*log(temp) - a + 1
			}else{
				#r <- A^(-2*(a-1)) * exp(- b*A^(-2) / 2)
			logr <- -2*a + 2 - b*A^2/2
			}
		}		
	}
	i <- 1
	j <- 0	
	if( source == 'gamma') {
	  repeat{
			x <- rgamma(1,a,b)
			j <- j + 1
			if( x > A^(-2)){
				out[i] <- x
				i <- i + 1
			}
			if( i > n) break
		}	
	}else{
		repeat{
		  ## assume source = exp
		  x <- A^(-2) + rexp(1,lambda)
		  j <- j + 1
		  lacc <- -logr + (a-1)*log(x) - (b -lambda)*x
		  u <- runif(1)
		  if( log(u) <= lacc) {
		  	out[i] <- x
		  	i <- i + 1
		  }
			if( i > n) break	
		}
	}
	# print(n / j)
	out	
}

gentnorm <- function(b,q)
{
 ## Algorithm due to Geweke 91
 ## sample from N_c(b,q) I( x> 0)
 astar <- (-b/q) / (sqrt(1/q))
 i <- 0
 if( astar > .45 ) {
 	r <- exp( - (astar^2) / 2 )
 	repeat{
 		## rejection envelop a + exp(a^{-1})
 		xstar <- astar + rexp(1,1/astar)
 		i <- i + 1
 		u <- rexp(1)
 		lacc <- - xstar^2 / 2 + xstar / astar - 1 - log(r)
 		if( u > - lacc) break
 	}
 } else {
 	repeat{
 		## rejection envelop N(0,1)
 		xstar <- rnorm(1)
 		i <- i + 1
 		if( xstar > astar ) break
 	}
 }
 x <- b/q + sqrt(1/q) * xstar
 cat('gentnorm debug:',i)
 x
}
isload <- function(pkg)
{
  checkpkg <- grep(pkg,search())
  length(checkpkg) > 0
}
