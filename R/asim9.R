
# Resampling site by site with column species richness fixed (subroutine isc)
# Equivalent to algorithm ISR in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by sequentially (column by column)
# assigning individuals randomly to each column with probabilities
# proportional to observed row abundance totals until the respective column
# total species richness is reached.
#
#                     Preserved?
#     Total abundance:
#   Total occurrences:    √
#         Occurrences:
#  Species abundances:
# Species occurrences:
#     Site abundances:
#     Site richnesses:    √
#
# speciesData is a species-by-site abundance matrix

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
