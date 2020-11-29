\name{asim6}
\title{
Resampling with total species richness fixed
}
\description{
Randomizes an abundance matrix by randomly assigning individuals to matrix
cells with probabilities proportional to observed row and column abundance
totals until total species richness is reached.
}
\usage{
asim6(speciesData)
}
\arguments{
  \item{speciesData}{
species-by-site abundance matrix
}
}
\details{
Resamples total abundance conditional on marginal sums and constrained to observed species
richness; equivalent to algorithm IR in Ulrich and Gotelli (2010), and to the Fortran
\code{subroutine ir} in Ulrich (2008).
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
May generate empty columns (sites).
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
amat <- rmat()
asim6(amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
