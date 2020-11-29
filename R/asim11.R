
# Resampling species by species with row abundances fixed (subroutine itr)
# Equivalent to algorithm ITR in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by sequentially (row after row) assigning
# individuals randomly to each row with probabilities proportional to
# observed column abundance totals until the respective total row abundance
# is reached.
#
#                     Preserved?
#     Total abundance:    √
#   Total occurrences:
#         Occurrences:
#  Species abundances:    √
# Species occurrences:
#     Site abundances:
#     Site richnesses:
#
# speciesData is a species-by-site abundance matrix

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
