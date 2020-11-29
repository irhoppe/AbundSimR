\name{asim4}
\title{
Resampling with occurrences and total abundance fixed
}
\description{
Randomizes an abundance matrix by assigning individuals to nonempty
cells, with probabilities proportional to observed row and column abundance
totals. Individuals are sequentially assigned to the matrix in this way until
the total (original) matrix abundance is reached.
}
\usage{
asim4(speciesData)
}
\arguments{
  \item{speciesData}{
species-by-site abundance matrix
}
}
\details{
Resamples total abundance conditional on marginal sums; equivalent to algorithm
OS in Ulrich and Gotelli (2010), and to the Fortran \code{subroutine oa} in Ulrich (2008).

Preserves total matrix abundance and occurrences, but allows both species and site abundances
(marginal sums) to vary in proportion to their original values.
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
asim4(amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
