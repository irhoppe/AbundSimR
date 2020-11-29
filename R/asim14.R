
# Resampling with row/column species richness and abundances fixed (subroutine iff)
# Equivalent to algorithm IF in Ulrich and Gotelli (2010).
#
# Two-step algorithm for randomizing an abundance matrix.
#
#                     Preserved?
#     Total abundance:    ?
#   Total occurrences:    ?
#         Occurrences:    ?
#  Species abundances:    ?
# Species occurrences:    ?
#     Site abundances:    ?
#     Site richnesses:    ?
#
# speciesData is a species-by-site abundance matrix

asim14 <- function(speciesData){

  stop("asim14 is untested/incomplete")

  occData <- speciesData > 0
  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ospp <- rowSums(occData)
  Ssit <- colSums(occData)
  Ntot <- sum(speciesData)
  Otot <- sum(occData)

  speciesData <- swap(occData + 0)



}
