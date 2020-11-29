
# Resampling with total species richness fixed (subroutine ir)
# Equivalent to algorithm IR in Ulrich and Gotelli (2010).
#
# Randomizes an abundance matrix by randomly assigning individuals to matrix
# cells with probabilities proportional to observed row and column abundance
# totals until total species richness is reached.
#
#                     Preserved?
#     Total abundance:
#   Total occurrences:
#         Occurrences:
#  Species abundances:
# Species occurrences:
#     Site abundances:
#     Site richnesses:
#
# CAUTION: May generate empty columns (sites).
#
# speciesData is a species-by-site abundance matrix

asim6 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Stot <- sum(Nspp>0)

  speciesData <- matrix(0, spp, sit)

  sspp <- vector("numeric", spp)
  s <- 0
  while( s < Stot ){
    i <- sample(spp, 1, prob=Nspp)
    j <- sample(sit, 1, prob=Nsit)
    speciesData[i,j] <- speciesData[i,j] + 1
    if( sspp[i]==0 ){
      sspp[i] <- 1
      s <- s + 1
      if( s > Stot ) break
    }
  }

  return(speciesData)

}
