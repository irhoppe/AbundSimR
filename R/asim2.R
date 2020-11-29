
asim2 <- function(speciesData){

  spp <- nrow(speciesData)
  sit <- ncol(speciesData)
  Ospp <- rowSums(speciesData>0)

  for( i in 1:spp ){
    if( Ospp[i] <= 1 ) next
    swaps <- 100*sit + rbinom(1,1,0.5)
    for( s in 1:swaps ){
      j <- sample(sit, 2, prob=speciesData[i,]>0)
      a <- speciesData[i,j[1]]
      speciesData[i,j[1]] <- speciesData[i,j[2]]
      speciesData[i,j[2]] <- a
    }
  }

  return(speciesData)

}
