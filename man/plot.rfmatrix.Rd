\name{plot.rfmatrix}
\alias{plot.rfmatrix}

\title{Plot recombination fractions or LOD scores for a single marker}

\description{
  Plot a slice (corresponding to a single marker) through the pairwise
  recombination fractions or LOD scores calculated by
  \code{\link{est.rf}} and extracted with \code{\link{pull.rf}}.
}

\usage{
\method{plot}{rfmatrix}(x, marker, \dots)
}

\arguments{
 \item{x}{An object of class \code{rfmatrix}, as output by \code{\link{pull.rf}}.}
 \item{marker}{A single marker name, as a character string.}
 \item{\dots}{Optional arguments passed to \code{\link{plot.scanone}}.}
}

\value{
An object of class \code{"scanone"} (as output by \code{\link{scanone}},
and which may be summarized by \code{\link{summary.scanone}} or plotted
with \code{\link{plot.scanone}}), containing the estimated recombination
fractions or LOD scores for the input marker against all others.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=5)}
fake.f2 <- est.rf(fake.f2)
marker <- markernames(fake.f2, chr=5)[6]
lod <- pull.rf(fake.f2, "lod")
plot(lod, marker, bandcol="gray70")
}

\seealso{ \code{\link{pull.rf}}, \code{\link{est.rf}}, \code{\link{plotRF}} }


\keyword{hplot}
