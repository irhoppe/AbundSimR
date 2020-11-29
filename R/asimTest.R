
#' Tool for testing randomization algorithms
#'
#' Test for preservation of various matrix metrics by abundance matrix randomization algorithms
#' (null models)
#'
#' @details
#'
#'
#' @param func [asim*()] function to test (passed as symbol)
#' @param tmat test matrix; default option generates a random matrix using [rmat()]
#' @param r number of randomizations to perform
#'
#' @seealso
#' [rmat()], [summarize.mat()]
#'
#' @export
#'
#' @examples
#' asimTest(asim10, r=100)
asimTest <- function( func, tmat=NULL, r=10, ... ){

  if( !is.null(tmat) ) tmat <- rmat(...)

  cat("\nPerforming randomizations...\n\n")
  nulls <- replicate(r, func(tmat))

  tmatsum <- summarize.mat(tmat)
  nullsums <- do.call( c, apply(nulls, 3, summarize.mat) )

  check <- data.frame( stat=names(tmatsum), preserved=NA_integer_ )
  checks <- length(tmatsum)

  for( i in 1:checks ){
    inds <- seq( from=i, to=r*checks, by=checks )
    nullvar <- nullsums[inds]
    tmatvar <- tmatsum[i][[1]]
    result <- suppressWarnings(sum(sapply( nullvar, function(.v){
      all(.v==tmatvar)
    })))
    check$preserved[i] <- result
  }

  cat("TEST MATRIX:\n")
  cat(sprintf("  Number of species: %6d\n", nrow(tmat)))
  cat(sprintf("    Number of sites: %6d\n", ncol(tmat)))
  cat(sprintf("  Total occurrences: %6d\n", tmatsum$Otot))
  cat(sprintf("    Total abundance: %6d\n\n", tmatsum$Ntot))

  cat("RANDOMIZATIONS:\n")
  cat(sprintf("         Iterations: %6d\n", r))
  cat(sprintf("          Algorithm: %6s\n\n", as.character(substitute(func))))

  cat("RESULTS:                #/R\n")
  cat(sprintf("    Total abundance: % 4d/%d\n", check$preserved[match("Ntot",check$stat)], r))
  cat(sprintf("  Total occurrences: % 4d/%d\n", check$preserved[match("Otot",check$stat)], r))
  cat(sprintf("        Occurrences: % 4d/%d\n", check$preserved[match("Occ",check$stat)], r))
  cat(sprintf(" Species abundances: % 4d/%d\n", check$preserved[match("Nspp",check$stat)], r))
  cat(sprintf("Species occurrences: % 4d/%d\n", check$preserved[match("Ospp",check$stat)], r))
  cat(sprintf("    Site abundances: % 4d/%d\n", check$preserved[match("Nsit",check$stat)], r))
  cat(sprintf("   Site occurrences: % 4d/%d\n", check$preserved[match("Osit",check$stat)], r))

}
