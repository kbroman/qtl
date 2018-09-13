\name{badorder}
\alias{badorder}
\docType{data}

\title{An intercross with misplaced markers}

\description{
  Simulated data for an intercross with some markers out of order.
}

\usage{data(badorder)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details.
}

\details{
  There are 250 F2 individuals typed at a total of 36 markers on four
  chromosomes.  The data were simulated with QTLs at the center of
  chromosomes 1 and 3.

  The order of several markers on chromosome 1 is incorrect.  Markers on
  chromosomes 2 and 3 are switched.
}

\examples{
data(badorder)
\dontshow{badorder <- subset(badorder,chr=1)
}
# estimate recombination fractions
badorder <- est.rf(badorder)
plotRF(badorder)

# re-estimate map
newmap <- est.map(badorder)
plotMap(badorder, newmap)

# assess marker order on chr 1
rip3 <- ripple(badorder, chr=1, window=3)
summary(rip3)
}

\seealso{ \code{\link{est.rf}}, \code{\link{ripple}},
  \code{\link{est.map}}, \code{\link{sim.cross}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{datasets}
