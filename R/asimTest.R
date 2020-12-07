
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
#' @param ... additional arguments to be passed to [rmat()]
#'
#' @seealso
#' [rmat()], [summarize.mat()]
asimTest <- function( func, tmat=NULL, rep=100, ... ){

  if( is.null(tmat) ) tmat <- rmat(...)

  cat("\nPerforming randomizations...\n\n")
  t0 <- Sys.time()
  nulls <- replicate(rep, func(tmat))
  t1 <- Sys.time()
  t <- t1-t0

  tmatsum <- summarize.mat(tmat)
  nullsums <- apply(nulls, 3, summarize.mat)

  Occ <- vector( "list", rep )
  Nspp <- Ospp <- matrix( 0, rep, nrow(tmat) )
  Nsit <- Osit <- matrix( 0, rep, ncol(tmat) )
  for( r in 1:rep ){
    Occ[[r]] <- nullsums[[r]]$Occ
    Nspp[r,] <- nullsums[[r]]$Nspp
    Ospp[r,] <- nullsums[[r]]$Ospp
    Nsit[r,] <- nullsums[[r]]$Nsit
    Osit[r,] <- nullsums[[r]]$Osit
  }
  nNspp <- colMeans(Nspp)/max(c(colMeans(Nspp),tmatsum$Nspp))
  nNsit <- colMeans(Nsit)/max(c(colMeans(Nsit),tmatsum$Nsit))
  nOspp <- colMeans(Ospp)/max(c(colMeans(Ospp),tmatsum$Ospp))
  nOsit <- colMeans(Osit)/max(c(colMeans(Osit),tmatsum$Osit))
  tNspp <- tmatsum$Nspp /max(c(colMeans(Nspp),tmatsum$Nspp))
  tNsit <- tmatsum$Nsit /max(c(colMeans(Nsit),tmatsum$Nsit))
  tOspp <- tmatsum$Ospp /max(c(colMeans(Ospp),tmatsum$Ospp))
  tOsit <- tmatsum$Osit /max(c(colMeans(Osit),tmatsum$Osit))

  par( mfrow=c(2,2), mar=c(4,3,2,1), mgp=c(1,0,0), family="serif", xaxt="n", yaxt="n" )
  plot( tNspp, xlab="",        ylab="Abundance"  , col="#8080FF" ); lines(nNspp)
  plot( tNsit, xlab="",        ylab=""           , col="#8080FF" ); lines(nNsit)
  plot( tOspp, xlab="Species", ylab="Occurrences", col="#8080FF" ); lines(nOspp)
  plot( tOsit, xlab="Site",    ylab=""           , col="#8080FF" ); lines(nOsit)
  mtext( as.character(substitute(func)), outer=TRUE )
  invisible(readline(prompt="Press [enter] to continue"))
  plot( sort(tNspp), xlab="",        ylab="Abundance"  , col="#8080FF" ); lines(nNspp[order(tNspp)])
  plot( sort(tNsit), xlab="",        ylab=""           , col="#8080FF" ); lines(nNsit[order(tNsit)])
  plot( sort(tOspp), xlab="Species", ylab="Occurrences", col="#8080FF" ); lines(nOspp[order(tOspp)])
  plot( sort(tOsit), xlab="Site",    ylab=""           , col="#8080FF" ); lines(nOsit[order(tOsit)])

  nullsums <- do.call( c, nullsums )

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
  cat(sprintf("          Algorithm: %6s\n", as.character(substitute(func))))
  cat(sprintf("       Elapsed time: %.2f %s\n\n", t, units(t)))

  cat("RESULTS:                #/R\n")
  cat(sprintf("    Total abundance: % 4d/%d\n", check$preserved[match("Ntot",check$stat)], rep))
  cat(sprintf("  Total occurrences: % 4d/%d\n", check$preserved[match("Otot",check$stat)], rep))
  cat(sprintf("        Occurrences: % 4d/%d\n", check$preserved[match( "Occ",check$stat)], rep))
  cat(sprintf(" Species abundances: % 4d/%d\n", check$preserved[match("Nspp",check$stat)], rep))
  cat(sprintf("Species occurrences: % 4d/%d\n", check$preserved[match("Ospp",check$stat)], rep))
  cat(sprintf("    Site abundances: % 4d/%d\n", check$preserved[match("Nsit",check$stat)], rep))
  cat(sprintf("   Site occurrences: % 4d/%d\n", check$preserved[match("Osit",check$stat)], rep))

}
