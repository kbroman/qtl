\name{pull.rf}
\alias{pull.rf}

\title{Pull out recombination fractions or LOD scores from a cross object}

\description{
  Pull out either the pairwise recombination fractions or the LOD
  scores, as calculated by \code{\link{est.rf}}, from a cross object.
}

\usage{pull.rf(cross, what=c("rf", "lod"), chr)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{what}{Indicates whether to pull out a matrix of estimated
 recombination fractions or a matrix of LOD scores.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
}

\value{
An object of class \code{"rfmatrix"}, which is a matrix of either
estimated recombination fractions between all marker pairs or of LOD
scores (for the test of rf=1/2) for all marker pairs.

The genetic map is included as an attribute.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=5)}
fake.f2 <- est.rf(fake.f2)
rf <- pull.rf(fake.f2)
lod <- pull.rf(fake.f2, "lod")
plot(rf[1,], lod[1,], xlab="rec frac", ylab="LOD score")
marker <- markernames(fake.f2, chr=5)[6]
par(mfrow=c(2,1))
plot(rf, marker, bandcol="gray70")
plot(lod, marker, bandcol="gray70")
}

\seealso{ \code{\link{est.rf}}, \code{\link{plot.rfmatrix}}, \code{\link{plotRF}} }


\keyword{utilities}
