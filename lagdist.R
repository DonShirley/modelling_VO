lagdist <- function(x,maxlag=round(dim(x)[1]/4)) {
  #----------------------------------------------
  # Estimation of mean square displacement
  # (Collins & DeLuca, 1993, Exp Brain Res)
  # R-Version 1.0, 28-OCT-2014 
  #----------------------------------------------
  x1 <- x
  x2 <- x
  r <- rep(0,maxlag)
  for ( l in 1:maxlag )  {
    x1 <- x1[-1,]
    x2 <- x2[-dim(x2)[1],]
    d <- x1-x2
    r[l] <- mean( d[,1]^2 + d[,2]^2 )
  }
  return(r)
}