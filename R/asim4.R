
# Resampling with occurrences and total abundance fixed (subroutine oa)
# Equivalent to algorithm OS in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by randomly assigning individuals to nonempty
# cells, with probabilities proportional to observed row and column abundance
# totals. Individuals are sequentially assigned to the matrix in this way until,
# the total (original) matrix abundance is reached.
#
#                     Preserved?
#     Total abundance:    √
#   Total occurrences:    √
#         Occurrences:    √
#  Species abundances:
# Species occurrences:    √
#     Site abundances:
#     Site richnesses:    √
#
# speciesData is a species-by-site abundance matrix

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
