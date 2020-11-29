\name{asim8}
\title{
Resampling species by species with row site occurrences fixed
}
\description{
Randomizes an abundance matrix by sequentially (row by row) assigning
individuals to each row with probabilities proportional to
observed column abundance totals until the respective number of row
occurrences is reached.
}
\usage{
asim8(speciesData)
}
\arguments{
  \item{speciesData}{
species-by-site abundance matrix
}
}
\details{
Resamples species abundances conditional on column sums and constrained to observed species
richness; equivalent to algorithm ISR in Ulrich and Gotelli (2010), and to the Fortran
\code{subroutine isr} in Ulrich (2008).

Preserves total occurrences and species occurrences, but allows abundances and site richnesses
to vary.
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
asim8(amat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
