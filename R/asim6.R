
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
