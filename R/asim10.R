
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
