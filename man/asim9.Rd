\name{asim9}
\title{
Resampling site by site with column species richness fixed
}
\description{
Randomizes an abundance matrix by sequentially (column by column)
assigning individuals to each column with probabilities proportional
to observed row abundance totals until the respective column total
species richness is reached.
}
\usage{
asim9(speciesData)
}
\arguments{
  \item{speciesData}{
species-by-site abundance matrix
}
}
\details{
Resamples species abundances conditional on row sums and constrained to observed
site occurrences; equivalent to algorithm ISC in Ulrich and Gotelli (2010),
and to the Fortran \code{subroutine isc} in Ulrich (2008).

Preserves total occurrences and site richnesses, but allows abundances and species
occurrences and to vary.
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
asim9(amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
