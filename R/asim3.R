
# Abundance swap within columns (subroutine pc)
# Equivalent to algorithm PC in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by reshuffling populations equiprobably
# among the nonempty cells of each column.
#
#                     Preserved?
#     Total abundance:    √
#   Total occurrences:    √
#         Occurrences:    √
#  Species abundances:
# Species occurrences:    √
#     Site abundances:    √
#     Site richnesses:    √
#
# speciesData is a species-by-site abundance matrix

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
