\name{asim2}
\title{
Abundance swap within rows
}
\description{
Randomizes an abundance matrix by reshuffling populations equiprobably among the nonempty cells of each row.
}
\usage{
asim2(speciesData)
}
\arguments{
  \item{speciesData}{
species-by-site abundance matrix
}
}
\details{
Performs a series of abundance swaps within rows (species); equivalent to algorithm PR
in Ulrich and Gotelli (2010), and to the Fortran \code{subroutine pr} in Ulrich (2008).

Preserves total matrix abundance and occurrences as well as species abundances (row totals),
but allows site abundances (column totals) to vary equiprobably.
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
asim2(amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
