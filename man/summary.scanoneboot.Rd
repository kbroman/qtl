\name{summary.scanoneboot}
\alias{summary.scanoneboot}

\title{Bootstrap confidence interval for QTL location}

\description{
  Calculates a bootstrap confidence interval for QTL location, using the
  bootstrap results from \code{\link{scanoneboot}}.
}

\usage{
\method{summary}{scanoneboot}(object, prob=0.95, expandtomarkers=FALSE, \dots)
}

\arguments{
 \item{object}{Output from \code{\link{scanoneboot}}.}
 \item{prob}{Desired coverage.}
 \item{expandtomarkers}{If TRUE, the interval is expanded to the nearest
   flanking markers.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{scanone}, indicating the
  position with the maximum LOD, and indicating endpoints
  for the estimated bootstrap confidence interval.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
\dontrun{data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2, step=1, err=0.001)
bootoutput <- scanoneboot(fake.f2, chr=13, method="hk")

summary(bootoutput)}
}

\seealso{ \code{\link{scanoneboot}}, \code{\link{plot.scanoneboot}},
  \code{\link{lodint}}, \code{\link{bayesint}} }

\keyword{utilities}
