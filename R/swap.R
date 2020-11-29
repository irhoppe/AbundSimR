swap <- function(occData){

  spp <- nrow(occData)
  sit <- ncol(occData)
  spsi <- length(occData)

  swaps <- 100*spsi + rbinom(1,1,0.5)
  s <- 0

  while( s < swaps ){
    print(s)
    i <- sample(spp, 2)
    j <- sample(sit, 2)
    if( (occData[i[1],j[1]] != 0 & occData[i[2],j[2]] !=0) | (occData[i[2],j[1]] != 0 & occData[i[1],j[2]] != 0) ){
      if( (occData[i[1],j[1]] == 0 & occData[i[2],j[2]] == 0) | (occData[i[2],j[1]] == 0 & occData[i[1],j[2]] == 0) ){
        s <- s + 1
        a <- occData[i[1],j[1]]
        occData[i[1],j[1]] <- occData[i[1],j[2]]
        occData[i[1],j[2]] <- a
        a <- occData[i[2],j[1]]
        occData[i[2],j[1]] <- occData[i[2],j[2]]
        occData[i[2],j[2]] <- a
      }
    }
  }

  return(occData)

}
