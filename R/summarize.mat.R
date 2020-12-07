summarize.mat <- function(.mat){
  .bmat <- .mat>0         # create a binary occurrence matrix
  list(
    Ntot=sum(.mat),       # total abundance
    Occ=which(.bmat),     # indices of occurrence
    Otot=sum(.bmat),      # total occurrences
    Nspp=rowSums(.mat),   # species abundances
    Ospp=rowSums(.bmat),  # species occurrences
    Nsit=colSums(.mat),   # site abundances
    Osit=colSums(.bmat)   # site occurrences
  )
}
