
# Resampling site by site with column abundances fixed (subroutine itc)
# Equivalent to algorithm ITC in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by sequentially (column after column)
# assigning individuals randomly to each column with probabilities
# proportional to observed row abundance totals until the respective total
# column abundance is reached.
#
#                     Preserved?
#     Total abundance:    √
#   Total occurrences:
#         Occurrences:
#  Species abundances:
# Species occurrences:
#     Site abundances:    √
#     Site richnesses:
#
# speciesData is a species-by-site abundance matrix

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
