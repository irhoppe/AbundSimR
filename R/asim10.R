
# Resampling with row/column abundances fixed (subroutine it)
# Equivalent to algorithm IT in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by randomly assigning individuals to matrix
# cells with probabilities proportional to observed row and column abundance
# totals until, for each row and column, total abundances are reached. Does not
# preserve species occurrence/non-occurrence.
#
#                     Preserved?
#     Total abundance:    √
#   Total occurrences:
#         Occurrences:
#  Species abundances:    √
# Species occurrences:
#     Site abundances:    √
#     Site richnesses:
#
# speciesData is a species-by-site abundance matrix

asim10 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ntot <- sum(speciesData)

  col <- vector("numeric", sit)
  row <- vector("numeric", spp)
  N <- 0
  speciesData <- matrix(0, spp, sit)

  while( N < Ntot ){
    i <- sample(1:spp, 1, prob=Nspp)
    while( row[i] >= Nspp[i] ) i <- sample(spp, 1, prob=Nspp)
    j <- sample(1:sit, 1, prob=Nsit)
    while( col[j] >= Nsit[j] ) j <- sample(sit, 1, prob=Nsit)
    row[i] <- row[i] + 1
    col[j] <- col[j] + 1
    speciesData[i,j] <- speciesData[i,j] + 1
    N <- N + 1
  }

  return(speciesData)

}
