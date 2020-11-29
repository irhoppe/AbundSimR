
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
