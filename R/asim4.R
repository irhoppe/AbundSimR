
asim4 <- function(speciesData){

  occData <- speciesData > 0
  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Nspp <- rowSums(speciesData)
  Nsit <- colSums(speciesData)
  Ntot <- sum(speciesData)
  Otot <- sum(occData)
  spsi <- length(speciesData)

  speciesData <- occData + 0

  while( Otot < Ntot ){
    ij <- sample(spsi, 1, prob=occData)
    speciesData[ij] <- speciesData[ij] + 1
    Otot <- Otot + 1
  }

  return(speciesData)

}
