\name{bayesint}
\alias{bayesint}

\title{Bayesian credible interval}

\description{
  Calculate an approximate Bayesian credible interval for a particular
  chromosome, using output from \code{\link{scanone}}.
}

\usage{bayesint(results, chr, qtl.index, prob=0.95, lodcolumn=1, expandtomarkers=FALSE)}

\arguments{
 \item{results}{Output from \code{\link{scanone}}, or a qtl object
 as output from \code{\link{refineqtl}}.}
 \item{chr}{A chromosome ID (if input \code{results} are from
 \code{\link{scanone}} (should have length 1).} 
 \item{qtl.index}{Numeric index for a QTL (if input \code{results} are
 from \code{\link{refineqtl}} (should have length 1).} 
 \item{prob}{Probability coverage of the interval.}
 \item{lodcolumn}{An integer indicating which
   of the LOD score columns should be considered (if input
 \code{results} are from \code{\link{scanone}}).}
 \item{expandtomarkers}{If TRUE, the interval is expanded to the nearest
   flanking markers.}
}

\value{
  An object of class \code{scanone} indicating the
  estimated QTL position and the approximate endpoints
  for the Bayesian credible interval.
}

\details{
  We take \eqn{10^{LOD}}{10^LOD}, rescale it to have area 1, and then
  calculate the connected interval with density above some threshold
  and having coverage matching the target probability.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,4))}
hyper <- calc.genoprob(hyper, step=0.5)
out <- scanone(hyper, method="hk")
bayesint(out, chr=1)
bayesint(out, chr=4)
bayesint(out, chr=4, prob=0.99)
bayesint(out, chr=4, expandtomarkers=TRUE)
}

\seealso{ \code{\link{scanone}}, \code{\link{lodint}} }

\keyword{utilities}
