
# Resampling species by species with row species richness fixed (subroutine isr)
# Equivalent to algorithm ISR in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by sequentially (row by row) assigning
# individuals randomly to each row with probabilities proportional to
# observed column abundance totals until the respective number of row
# occurrences is reached.
#
#                     Preserved?
#     Total abundance:
#   Total occurrences:    √
#         Occurrences:
#  Species abundances:
# Species occurrences:    √
#     Site abundances:
#     Site richnesses:
#
# speciesData is a species-by-site abundance matrix

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
