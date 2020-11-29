
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
