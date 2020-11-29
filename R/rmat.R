rmat <- function( nr=NULL, nc=NULL, p=NULL, N=NULL ){

  if( is.null(nr) ) nr <- round(runif(1,5,25))
  if( is.null(nc) ) nc <- round(runif(1,5,25))
  if( is.null(p) ) p <- runif(1)
  if( is.null(N) ) N <- 12
  S <- nr*nc
  cells <- rbinom(S,1,p)*rpois(S,N)

  return(matrix(cells,nr,nc))

}
