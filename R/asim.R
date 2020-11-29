
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
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim1(amat)         # performs equiprobable abundance swap algorithm
asim1 <- function(speciesData){

  occData <- speciesData > 0
  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  spsi <- length(speciesData)

  swaps <- 100*spsi + rbinom(1,1,0.5)

  for( s in 1:swaps ){
    ij <- sample(spsi, 2, prob=occData)
    a <- speciesData[ij[1]]
    speciesData[ij[1]] <- speciesData[ij[2]]
    speciesData[ij[2]] <- a
  }

  return(speciesData)

}

#' Abundance swap within rows
#'
#' Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells
#' of each row.
#'
#' @details
#' Performs a series of abundance swaps within rows (species); equivalent to algorithm PR in Ulrich
#' and Gotelli (2010), and to the Fortran `subroutine pr` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences as well as species abundances (row totals), but
#' allows site abundances (column totals) to vary equiprobably.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim2(amat)         # performs equiprobable row abundance swap algorithm
asim2 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Ospp <- rowSums(speciesData>0)

  for( i in 1:spp ){
    if( Ospp[i] <= 1 ) next
    swaps <- 100*sit + rbinom(1,1,0.5)
    for( s in 1:swaps ){
      j <- sample(sit, 2, prob=speciesData[i,]>0)
      a <- speciesData[i,j[1]]
      speciesData[i,j[1]] <- speciesData[i,j[2]]
      speciesData[i,j[2]] <- a
    }
  }

  return(speciesData)

}

#' Abundance swap within columns
#'
#' Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells
#' of each column.
#'
#' @details
#' Performs a series of abundance swaps within columns (sites); equivalent to algorithm PC in
#' Ulrich and Gotelli (2010), and to the Fortran `subroutine pc` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences as well as site abundances (column totals), but
#' allows species abundances (row totals) to vary equiprobably.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim3(amat)         # performs equiprobable column abundance swap algorithm
asim3 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Ssit <- colSums(speciesData>0)

  for( j in 1:sit ){
    if( Ssit[j] <= 1 ) next
    swaps <- 100*spp + rbinom(1,1,0.5)
    for( s in 1:swaps ){
      i <- sample(spp, 2, prob=speciesData[,j]>0)
      a <- speciesData[i[1],j]
      speciesData[i[1],j] <- speciesData[i[2],j]
      speciesData[i[2],j] <- a
    }
  }

  return(speciesData)

}

#' Resampling with occurrences and total abundance fixed
#'
#' Randomizes an abundance matrix by assigning individuals to nonempty cells, with probabilities
#' proportional to observed row and column abundance totals. Individuals are sequentially assigned
#' to the matrix in this way until the total (original) matrix abundance is reached.
#'
#' @details
#' Resamples total abundance conditional on marginal sums; equivalent to algorithm OS in Ulrich and
#' Gotelli (2010), and to the Fortran `subroutine oa` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences, but allows both species and site abundances
#' (marginal sums) to vary in proportion to their original values.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim4(amat)         # performs proportional abundance resampling algorithm
asim4 <- function(speciesData){

  occData <- speciesData > 0
  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ntot <- sum(speciesData)
  Otot <- sum(occData)
  spsi <- length(speciesData)

  speciesData <- occData + 0

  while( Otot < Ntot ){
    ij <- sample(spsi, 1, prob=occData)
    speciesData[ij] <- speciesData[ij] + 1
    Otot <- Otot + 1
  }

  return(speciesData)

}

#' Resampling with occurrences and row/column abundances fixed
#'
#' Randomizes an abundance matrix by assigning individuals to nonempty cells, with probabilities
#' proportional to observed row and column abundance totals. Individuals are sequentially assigned
#' to the matrix in this way until, for each row and column, total abundances are reached.
#'
#' @details
#' Resamples total abundance conditional on and constrained by marginal sums; equivalent to algorithm
#' OF in Ulrich and Gotelli (2010), and to the Fortran `subroutine oaf` in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim5(amat)         # performs proportional abundance resampling algorithm
asim5 <- function(speciesData){

  stop("asim5 is untested/incomplete")

  occData <- speciesData > 0
  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ospp <- rowSums(occData)
  Ssit <- colSums(occData)
  Ntot <- sum(speciesData)
  Otot <- sum(occData)

  speciesData <- occData + 0
  r <- 0

  while( Otot < Ntot ){
    i <- sample(spp, 1, prob=Nspp)
    j <- sample(sit, 1, prob=Nsit)
    if( speciesData[i,j] > 0 ){
      r <- r + 1
      if( r > 10*Ntot | Otot == (Ntot-1) ) stop("Trials exceeded maximum")   # <--- CONDITION NEEDS WORK
      if( !(Ospp[i] >= Nspp[i] | Ssit[j] >= Nsit[j]) ){
        Ospp[i] <- Ospp[i] + 1
        Ssit[j] <- Ssit[j] + 1
        speciesData[i,j] <- speciesData[i,j] + 1
        Otot <- Otot + 1
      }
    }
  }

  return(speciesData)

}

#' Resampling with total species richness fixed
#'
#' Randomizes an abundance matrix by randomly assigning individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until total species richness is reached.
#'
#' @details
#' Resamples total abundance conditional on marginal sums and constrained to observed species
#' richness; equivalent to algorithm IR in Ulrich and Gotelli (2010), and to the Fortran
#' `subroutine ir` in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim6(amat)         # performs proportional abundance resampling algorithm
asim6 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Stot <- sum(Nspp>0)

  speciesData <- matrix(0, spp, sit)

  sspp <- vector("numeric", spp)
  s <- 0
  while( s < Stot ){
    i <- sample(spp, 1, prob=Nspp)
    j <- sample(sit, 1, prob=Nsit)
    speciesData[i,j] <- speciesData[i,j] + 1
    if( sspp[i]==0 ){
      sspp[i] <- 1
      s <- s + 1
      if( s > Stot ) break
    }
  }

  return(speciesData)

}

#' Resampling with row/column species richness fixed
#'
#' Randomizes an abundance matrix by assigning individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until the total number of occurrences
#' is reached for each row and column.
#'
#' @details
#' Resamples total abundance conditional on marginal sums and constrained to observed species
#' occurrences and site richnesses; equivalent to algorithm IS in Ulrich and Gotelli (2010), and to
#' the Fortran `subroutine is` in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim7(amat)         # performs proportional abundance resampling algorithm
asim7 <- function(speciesData){

  stop("asim7 is untested/incomplete")

  occData <- speciesData > 0
  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ospp <- rowSums(occData)
  Ssit <- rowSums(occData)
  Ntot <- sum(speciesData)
  Otot <- sum(occData)

  speciesData <- matrix(0, spp, sit)


  ospp <- vector("numeric", spp)
  ssit <- vector("numeric", sit)
  o <- n <- 0
  while( o < Otot ){
    i <- sample(spp, 1, prob=Nspp)
    j <- sample(sit, 1, prob=Nsit)
    n <- n + 1
    if( n > (100*Ntot) ) break
    if( speciesData[i,j]==0 ){
      if( (ospp[i] > Ospp[i]) | (ssit[j] > Ssit[j]) ) next
      ospp[i] <- ospp[i] + 1
      ssit[j] <- ssit[j] + 1
      speciesData[i,j] <- speciesData[i,j] + 1
      o <- o + 1
    } else {
      speciesData[i,j] <- speciesData[i,j] + 1
    }
  }

  return(speciesData)

}

#' Resampling species by species with row site occurrences fixed
#'
#' Randomizes an abundance matrix by sequentially (row by row) assigning individuals to each row
#' with probabilities proportional to observed column abundance totals until the respective
#' number of row occurrences is reached.
#'
#' @details
#' Resamples species abundances conditional on column sums and constrained to observed species
#' richness; equivalent to algorithm ISR in Ulrich and Gotelli (2010), and to the Fortran
#' `subroutine isr` in Ulrich (2008).
#'
#' Preserves total occurrences and species occurrences, but allows abundances and site richnesses
#' to vary.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim8(amat)         # performs proportional abundance resampling algorithm
asim8 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ospp <- rowSums(speciesData>0)

  speciesData <- matrix(0, spp, sit)

  for( i in 1:spp ){
    o <- 0
    while( o < Ospp[i] ){
      j <- sample(sit, 1, prob=Nsit)
      if( speciesData[i,j]==0 ) o <- o + 1
      speciesData[i,j] <- speciesData[i,j] + 1
    }
  }

  return(speciesData)

}

#' Resampling site by site with column species richness fixed
#'
#' Randomizes an abundance matrix by sequentially (column by column) assigning individuals to each
#' column with probabilities proportional to observed row abundance totals until the respective
#' column total species richness is reached.
#'
#' @details
#' Resamples species abundances conditional on row sums and constrained to observed site occurrences;
#' equivalent to algorithm ISC in Ulrich and Gotelli (2010), and to the Fortran `subroutine isc` in
#' Ulrich (2008).
#'
#' Preserves total occurrences and site richnesses, but allows abundances and species occurrences to
#' vary.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim9(amat)         # performs proportional abundance resampling algorithm
asim9 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ssit <- colSums(speciesData>0)

  speciesData <- matrix(0, spp, sit)

  for( j in 1:sit ){
    s <- 0
    while( s < Ssit[j] ){
      i <- sample(spp, 1, prob=Nspp)
      if( speciesData[i,j]==0 ) s <- s + 1
      speciesData[i,j] <- speciesData[i,j] + 1
    }
  }

  return(speciesData)

}

#' Resampling with row/column abundances fixed
#'
#' Randomizes an abundance matrix by assigning individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until, for each row and column,
#' total abundances are reached.
#'
#' @details
#' Resamples species abundances conditional on and constrained to row and column sums and
#' constrained to observed species richness; equivalent to algorithm IT in Ulrich and Gotelli
#' (2010), and to the Fortran `subroutine itt` in Ulrich (2008).
#'
#' Preserves species and site abundances, but allows occurrences to vary.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim10(amat)        # performs proportional abundance resampling algorithm
asim10 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ntot <- sum(speciesData)

  col <- vector("numeric", sit)
  row <- vector("numeric", spp)
  N <- 0
  speciesData <- matrix(0, spp, sit)

  while( N < Ntot ){
    i <- sample(1:spp, 1, prob=Nspp)
    while( row[i] >= Nspp[i] ) i <- sample(spp, 1, prob=Nspp)
    j <- sample(1:sit, 1, prob=Nsit)
    while( col[j] >= Nsit[j] ) j <- sample(sit, 1, prob=Nsit)
    row[i] <- row[i] + 1
    col[j] <- col[j] + 1
    speciesData[i,j] <- speciesData[i,j] + 1
    N <- N + 1
  }

  return(speciesData)

}

#' Resampling species by species with row abundances fixed
#'
#' Randomizes an abundance matrix by sequentially (row after row) assigning individuals to each
#' row with probabilities proportional to observed column abundance totals until the respective
#' total row abundance is reached.
#'
#' @details
#' Resamples species abundances conditional on column sums and constrained to observed row abundance;
#' equivalent to algorithm ITR in Ulrich and Gotelli (2010), and to the Fortran `subroutine itr` in
#' Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim11(amat)        # performs proportional abundance resampling algorithm
asim11 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)

  speciesData <- matrix(0, spp, sit)

  for( i in 1:spp ){
    n <- 0
    while( n < Nspp[i] ){
      j <- sample(sit, 1, prob=Nsit)
      speciesData[i,j] <- speciesData[i,j] + 1
      n <- n + 1
    }
  }

  return(speciesData)

}

#' Resampling site by site with column abundances fixed
#'
#' Randomizes an abundance matrix by sequentially (column after column) assigning individuals to
#' each column with probabilities proportional to observed row abundance totals until the respective
#' total column abundance is reached.
#'
#' @details
#' Resamples species abundances conditional on row sums and constrained to observed column abundance;
#' equivalent to algorithm ITC in Ulrich and Gotelli (2010), and to the Fortran `subroutine itc` in
#' Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim12(amat)        # performs proportional abundance resampling algorithm
asim12 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)

  speciesData <- matrix(0, spp, sit)

  for( j in 1:sit ){
    n <- 0
    while( n < Nsit[j] ){
      i <- sample(spp, 1, prob=Nspp)
      speciesData[i,j] <- speciesData[i,j] + 1
      n <- n + 1
    }
  }

  return(speciesData)

}

#' Resampling with total abundance fixed
#'
#' Randomizes an abundance matrix by assigning all individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until the matrix-wide total abundance
#' is reached.
#'
#' @details
#' Resamples species abundances conditional on row and column sums and constrained to observed total
#' abundance; equivalent to algorithm IA in Ulrich and Gotelli (2010), and to the Fortran `subroutine ia`
#' in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty rows(species) or columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim13(amat)        # performs proportional abundance resampling algorithm
asim13 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ntot <- sum(speciesData)

  n <- 0

  speciesData <- matrix(0, spp, sit)

  while( n < Ntot ){
    i <- sample(spp, 1, prob=Nspp)
    j <- sample(sit, 1, prob=Nsit)
    speciesData[i,j] <- speciesData[i,j] + 1
    n <- n + 1
  }

  return(speciesData)

}

#' Resampling with row/column species richness and abundances fixed
#'
#' Two-step algorithm for randomizing an abundance matrix.
#'
#' @details
#' Equivalent to algorithm IF in Ulrich and Gotelli (2010), and to the Fortran `subroutine iff` in
#' Gotelli (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty rows(species) or columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20×12 abundance matrix
#' asim14(amat)        # performs proportional abundance resampling algorithm
asim14 <- function(speciesData){

  stop("asim14 is untested/incomplete")

  occData <- speciesData > 0
  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ospp <- rowSums(occData)
  Ssit <- colSums(occData)
  Ntot <- sum(speciesData)
  Otot <- sum(occData)

  speciesData <- swap(occData + 0)



}
