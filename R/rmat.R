
#' Random matrix generator
#'
#' Create a random abundance matrix, with optional controls for dimensions, probability of occurrence,
#' and mean cell abundance.
#'
#' @details
#'
#'
#' @param nr integer value for number of rows
#' @param nc integer value for number of columns
#' @param p proportion of cells to be non-zero
#' @param N integer value for mean value of non-zero cells
#'
#' @export
#'
#' @examples
#' nr <- 20
#' nc <- 12
#' p <- 0.35
#' N <- 20
#' amat <- rmat(nr,nc,p,N)   # creates a random 20 x 12 abundance matrix
#' dim(amat)                 # [1] 20 12 (nr nc)
#' sum(amat!=0)/length(amat) # converges to 0.35 (p)
#' mean(amat[amat!=0])       # converges to 20 (N)
rmat <- function( nr=NULL, nc=NULL, p=NULL, N=NULL ){

  if( is.null(nr) ) nr <- round(runif(1,5,25))
  if( is.null(nc) ) nc <- round(runif(1,5,25))
  if( is.null(p) ) p <- runif(1)
  if( is.null(N) ) N <- 12
  S <- nr*nc
  cells <- rbinom(S,1,p)*rpois(S,N)

  return(matrix(cells,nr,nc))

}
