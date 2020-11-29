\name{asim1}
\title{
Abundance swap
}
\description{
Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells
of the entire matrix.
}
\usage{
asim1(speciesData)
}
\arguments{
  \item{speciesData}{
species-by-site abundance matrix
}
}
\details{
Performs an abundance swap equivalent to algorithm PM in Ulrich and Gotelli (2010), and to
the Fortran \code{subroutine pm} in Ulrich (2008).

Preserves total matrix abundance and occurrences, but allows marginal sums (i.e. species and
site abundances) to vary equiprobably.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
amat <- rmat()
asim1(amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
