
#' Abundance swap
#'
#' Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells
#' of the entire matrix.
#'
#' @details
#' Performs an abundance swap equivalent to algorithm PM in Ulrich and Gotelli (2010), and to the
#' Fortran `subroutine pm` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences, but allows marginal sums (i.e. species and
#' site abundances) to vary equiprobably.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim1_f(amat)       # performs equiprobable abundance swap algorithm
asim1_f <- function(speciesData, rep=1L){

  stopifnot(is.matrix(speciesData), is.numeric(speciesData), is.numeric(rep))
  storage.mode(speciesData) <- "integer"

  shuffle <- .Fortran("pm", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), rep=as.integer(rep), PACKAGE="cooccur")
  return(shuffle$mat)

}
