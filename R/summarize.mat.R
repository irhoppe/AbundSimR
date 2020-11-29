summarize.mat <- function(.mat){
  .bmat <- .mat>0
  list(
    Ntot=sum(.mat),
    Occ=which(.bmat),
    Otot=sum(.bmat),
    Nspp=rowSums(.mat),
    Ospp=rowSums(.bmat),
    Nsit=colSums(.mat),
    Osit=colSums(.bmat)
  )
}
