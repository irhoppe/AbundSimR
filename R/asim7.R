
# Resampling with row/column species richness fixed (subroutine is)
# Equivalent to algorithm IS in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by randomly assigning individuals to matrix
# cells with probabilities proportional to observed row and column abundance
# totals until the total number of occurrences is reached for each row and
# column.
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
