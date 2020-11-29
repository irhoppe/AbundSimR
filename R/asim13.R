
# Resampling with total abundance fixed (subroutine ia (rand=='ia'))
# Equivalent to algorithm IA in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by randomly assigning all individuals to
# matrix cells with probabilities proportional to observed row and column
# abundance totals until the matrix-wide total abundance is reached.
# CAUTION: May generate empty rows (species) or columns (sites).
#
#                     Preserved?
#     Total abundance:    √
#   Total occurrences:
#         Occurrences:
#  Species abundances:
# Species occurrences:
#     Site abundances:
#     Site richnesses:
#
# speciesData is a species-by-site abundance matrix

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
