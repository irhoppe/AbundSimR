\name{asim10}
\title{
Resampling with row/column abundances fixed
}
\description{
Randomizes an abundance matrix by assigning individuals to matrix cells
with probabilities proportional to observed row and column abundance totals
until, for each row and column, total abundances are reached.
}
\usage{
asim10(speciesData)
}
\arguments{
  \item{speciesData}{
species-by-site abundance matrix
}
}
\details{
Resamples species abundances conditional on and constrained to row and column sums
and constrained to observed species richness; equivalent to algorithm IT in Ulrich
and Gotelli (2010), and to the Fortran \code{subroutine itt} in Ulrich (2008).

Preserves species and site abundances, but allows occurrences to vary.
}
\value{
Returns an abundance matrix of equal dimension to the input matrix.
}
\references{
Ulrich, W. 2008. Program Co-Occurrence. Version 1.

Ulrich, W. and Gotelli, N.J. 2010. Null model analysis of species associations using abundance data. Ecology 91:3384-3397.
}
\note{
%%  ~~further notes~~
}
\section{Caution }{
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
amat <- rmat()
asim10(amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
