\name{lodint}
\alias{lodint}

\title{LOD support interval}

\description{
  Calculate a LOD support interval for a particular chromosome,
  using output from scanone.
}

\usage{lodint(results, chr, qtl.index, drop=1.5, lodcolumn=1, expandtomarkers=FALSE)}

\arguments{
 \item{results}{Output from \code{\link{scanone}}, or a qtl object
 as output from \code{\link{refineqtl}}.}
 \item{chr}{A chromosome ID (if input \code{results} are from
 \code{\link{scanone}} (should have length 1).}
 \item{qtl.index}{Numeric index for a QTL (if input \code{results} are
 from \code{\link{refineqtl}} (should have length 1).}
 \item{drop}{LOD units to drop to form the interval.}
 \item{lodcolumn}{An integer indicating which
   of the LOD score columns should be considered (if input
 \code{results} are from \code{\link{scanone}}).}
 \item{expandtomarkers}{If TRUE, the interval is expanded to the nearest
   flanking markers.}
}

\value{
  An object of class \code{scanone} indicating the
  estimated QTL position and the approximate endpoints
  for the LOD support interval.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

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

\seealso{ \code{\link{scanone}}, \code{\link{bayesint}} }

\keyword{utilities}
