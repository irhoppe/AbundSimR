
# Abundance swap (subroutine pm)
# Equivalent to algorithm PM in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by reshuffling populations equiprobably
# among the nonempty cells of the entire matrix.
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
