
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
