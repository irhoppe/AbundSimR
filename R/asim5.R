
# Resampling with occurrences and row/column abundances fixed (subroutine of)
# Equivalent to algorithm OF in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by randomly assigning individuals to nonempty
# cells, with probabilities proportional to observed row and column abundance
# totals. Individuals are sequentially assigned to the matrix in this way until,
# for each row and column, total abundances are reached.
#
#                     Preserved?
#     Total abundance:    ?
#   Total occurrences:    ?
#         Occurrences:    ?
#  Species abundances:    ?
# Species occurrences:    ?
#     Site abundances:    ?
#     Site richnesses:    ?
#
# speciesData is a species-by-site abundance matrix

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
