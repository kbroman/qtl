\name{lodint}
\alias{lodint}

\title{LOD support interval}

\description{
  Calculate a LOD support interval for a particular chromosome,
  using output from scanone.
}

\usage{lodint(results, chr, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)}

\arguments{
 \item{results}{Output from \code{\link[qtl]{scanone}}.}
 \item{chr}{A chromosome ID.}
 \item{drop}{LOD units to drop to form the interval.}
 \item{lodcolumn}{An integer, or vector of 3 integers, indicating which
   of the LOD score columns should be plotted (generally this is 1).}
 \item{expandtomarkers}{If TRUE, the interval is expanded to the nearest
   flanking markers.}
}

\value{
  An object of class \code{scanone}, like the input, indicating the
  position with the maximum LOD, and indicating approximate endpoints
  for the LOD support interval.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,4))}
hyper <- calc.genoprob(hyper, step=0.5)
out <- scanone(hyper, method="hk")
lodint(out, chr=1)
lodint(out, chr=4)
lodint(out, chr=4, drop=2)
lodint(out, chr=4, expandtomarkers=TRUE)
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{bayesint}} }

\keyword{utilities}
