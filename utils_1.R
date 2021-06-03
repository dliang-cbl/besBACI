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