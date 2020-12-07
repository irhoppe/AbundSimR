
#' Abundance swap
#'
#' Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells
#' of the entire matrix.
#'
#' @details
#' Performs an abundance swap equivalent to algorithm PM in Ulrich and Gotelli (2010), and to the
#' Fortran `subroutine pm` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences, but allows marginal sums (i.e. species and
#' site abundances) to vary equiprobably.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim1_f(amat)       # performs equiprobable abundance swap algorithm
asim1_f <- function(speciesData, swaps=NULL){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  if( is.null(swaps) ){
    swaps <- 100*prod(dim(speciesData)) + rbinom(1,1,0.5)
  } else if( !is.numeric(swaps) ) stop("if provided, swaps must be a numeric value")

  shuffle <- .Fortran("pm", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), rep=as.integer(swaps[1]), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Abundance swap within columns
#'
#' Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells
#' of each column.
#'
#' @details
#' Performs a series of abundance swaps within columns (sites); equivalent to algorithm PC in
#' Ulrich and Gotelli (2010), and to the Fortran `subroutine pc` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences as well as site abundances (column totals), but
#' allows species abundances (row totals) to vary equiprobably.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim2_f(amat)       # performs equiprobable column abundance swap algorithm
asim2_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("pc", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Abundance swap within rows
#'
#' Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells
#' of each row.
#'
#' @details
#' Performs a series of abundance swaps within rows (species); equivalent to algorithm PR in Ulrich
#' and Gotelli (2010), and to the Fortran `subroutine pr` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences as well as species abundances (row totals), but
#' allows site abundances (column totals) to vary equiprobably.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim3_f(amat)       # performs equiprobable asrow abundance swap algorithm
asim3_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("pr", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling with occurrences and total abundance fixed
#'
#' Randomizes an abundance matrix by assigning individuals to nonempty cells, with probabilities
#' proportional to observed row and column abundance totals. Individuals are sequentially assigned
#' to the matrix in this way until the total (original) matrix abundance is reached.
#'
#' @details
#' Resamples total abundance conditional on marginal sums; equivalent to algorithm OS in Ulrich and
#' Gotelli (2010), and to the Fortran `subroutine oa` in Ulrich (2008).
#'
#' Preserves total matrix abundance and occurrences, but allows both species and site abundances
#' (marginal sums) to vary in proportion to their original values.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim4_f(amat)       # performs proportional abundance resampling algorithm
asim4_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("pr", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling with occurrences and row/column abundances fixed
#'
#' Randomizes an abundance matrix by assigning individuals to nonempty cells, with probabilities
#' proportional to observed row and column abundance totals. Individuals are sequentially assigned
#' to the matrix in this way until, for each row and column, total abundances are reached.
#'
#' @details
#' Resamples total abundance conditional on and constrained by marginal sums; equivalent to algorithm
#' OF in Ulrich and Gotelli (2010), and to the Fortran `subroutine of` in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim5_f(amat)       # performs proportional abundance resampling algorithm
asim5_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("of", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling with total species richness fixed
#'
#' Randomizes an abundance matrix by randomly assigning individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until total species richness is reached.
#'
#' @details
#' Resamples total abundance conditional on marginal sums and constrained to observed species
#' richness; equivalent to algorithm IR in Ulrich and Gotelli (2010), and to the Fortran
#' `subroutine ir` in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim6_f(amat)       # performs proportional abundance resampling algorithm
asim6_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("ir", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling with row/column species richness fixed
#'
#' Randomizes an abundance matrix by assigning individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until the total number of occurrences
#' is reached for each row and column.
#'
#' @details
#' Resamples total abundance conditional on marginal sums and constrained to observed species
#' occurrences and site richnesses; equivalent to algorithm IS in Ulrich and Gotelli (2010), and to
#' the Fortran `subroutine is` in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim7_f(amat)       # performs proportional abundance resampling algorithm
asim7_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("is", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling species by species with row site occurrences fixed
#'
#' Randomizes an abundance matrix by sequentially (row by row) assigning individuals to each row
#' with probabilities proportional to observed column abundance totals until the respective
#' number of row occurrences is reached.
#'
#' @details
#' Resamples species abundances conditional on column sums and constrained to observed species
#' richness; equivalent to algorithm ISR in Ulrich and Gotelli (2010), and to the Fortran
#' `subroutine isr` in Ulrich (2008).
#'
#' Preserves total occurrences and species occurrences, but allows abundances and site richnesses
#' to vary.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim8_f(amat)       # performs proportional abundance resampling algorithm
asim8_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("isr", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling site by site with column species richness fixed
#'
#' Randomizes an abundance matrix by sequentially (column by column) assigning individuals to each
#' column with probabilities proportional to observed row abundance totals until the respective
#' column total species richness is reached.
#'
#' @details
#' Resamples species abundances conditional on row sums and constrained to observed site occurrences;
#' equivalent to algorithm ISC in Ulrich and Gotelli (2010), and to the Fortran `subroutine isc` in
#' Ulrich (2008).
#'
#' Preserves total occurrences and site richnesses, but allows abundances and species occurrences to
#' vary.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim9_f(amat)       # performs proportional abundance resampling algorithm
asim9_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("isc", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling with row/column abundances fixed
#'
#' Randomizes an abundance matrix by assigning individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until, for each row and column,
#' total abundances are reached.
#'
#' @details
#' Resamples species abundances conditional on and constrained to row and column sums and
#' constrained to observed species richness; equivalent to algorithm IT in Ulrich and Gotelli
#' (2010), and to the Fortran `subroutine itt` in Ulrich (2008).
#'
#' Preserves species and site abundances, but allows occurrences to vary.
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim10_f(amat)      # performs proportional abundance resampling algorithm
asim10_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("itt", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling species by species with row abundances fixed
#'
#' Randomizes an abundance matrix by sequentially (row after row) assigning individuals to each
#' row with probabilities proportional to observed column abundance totals until the respective
#' total row abundance is reached.
#'
#' @details
#' Resamples species abundances conditional on column sums and constrained to observed row abundance;
#' equivalent to algorithm ITR in Ulrich and Gotelli (2010), and to the Fortran `subroutine itr` in
#' Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim11_f(amat)      # performs proportional abundance resampling algorithm
asim11_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("itr", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling site by site with column abundances fixed
#'
#' Randomizes an abundance matrix by sequentially (column after column) assigning individuals to
#' each column with probabilities proportional to observed row abundance totals until the respective
#' total column abundance is reached.
#'
#' @details
#' Resamples species abundances conditional on row sums and constrained to observed column abundance;
#' equivalent to algorithm ITC in Ulrich and Gotelli (2010), and to the Fortran `subroutine itc` in
#' Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim12_f(amat)      # performs proportional abundance resampling algorithm
asim12_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("itc", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling with total abundance fixed
#'
#' Randomizes an abundance matrix by assigning all individuals to matrix cells with probabilities
#' proportional to observed row and column abundance totals until the matrix-wide total abundance
#' is reached.
#'
#' @details
#' Resamples species abundances conditional on row and column sums and constrained to observed total
#' abundance; equivalent to algorithm IA in Ulrich and Gotelli (2010), and to the Fortran `subroutine ia`
#' in Ulrich (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty rows(species) or columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim13_f(amat)      # performs proportional abundance resampling algorithm
asim13_f <- function(speciesData){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  shuffle <- .Fortran("ia", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
  return(shuffle$mat)

}

#' Resampling with row/column species richness and abundances fixed
#'
#' Two-step algorithm for randomizing an abundance matrix.
#'
#' @details
#' Equivalent to algorithm IF in Ulrich and Gotelli (2010), and to the Fortran `subroutine iff` in
#' Gotelli (2008).
#'
#' @param speciesData species-by-site abundance matrix
#'
#' @return shuffled abundance matrix
#'
#' @section Caution:
#' May generate matrices with empty rows(species) or columns (sites).
#'
#' @family abundance matrix null models
#'
#' @references
#' Ulrich, W. 2008. Program Co-Occurrence. Version 1.
#'
#' Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. *Ecology* 91:3384-3397.
#'
#' @export
#'
#' @examples
#' amat <- rmat(20,12) # creates a random 20 x 12 abundance matrix
#' asim14_f(amat)      # performs proportional abundance resampling algorithm
asim14_f <- function(speciesData, swaps=NULL){

  if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
  if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"

  if( is.null(swaps) ){
    swaps <- 100*prod(dim(speciesData)) + rbinom(1,1,0.5)
  } else if( !is.numeric(swaps) ) stop("if provided, swaps must be a numeric value")

  shuffle <- .Fortran("iff", mat=speciesData, sp=as.integer(nrow(speciesData)),
                      si=as.integer(ncol(speciesData)), rep=as.integer(swaps[1]), PACKAGE="cooccur")
  return(shuffle$mat)

}

# asim15_f <- function(speciesData){
#
#   if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
#   if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"
#
#   shuffle <- .Fortran("o1", mat=speciesData, sp=as.integer(nrow(speciesData)),
#                       si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
#   return(shuffle$mat)
#
# }
#
# asim16_f <- function(speciesData){
#
#   if( !is.matrix(speciesData) | !is.numeric(speciesData ) ) stop("speciesData must be a numeric matrix")
#   if( !is.double(speciesData) )  storage.mode(speciesData) <- "double"
#
#   shuffle <- .Fortran("isa", mat=speciesData, sp=as.integer(nrow(speciesData)),
#                       si=as.integer(ncol(speciesData)), PACKAGE="cooccur")
#   return(shuffle$mat)
#
# }
