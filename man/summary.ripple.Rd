\name{summary.ripple}
\alias{summary.ripple}

\title{Print summary of ripple results}

\description{
  Print marker orders, from the output of the function \code{ripple},
  for which the log10 likelihood relative to the initial order is above
  a specified cutoff.
}

\usage{
\method{summary}{ripple}(object, lod.cutoff = -1, \dots)
}

\arguments{
 \item{object}{An object of class \code{ripple}, the output of the
   function \code{\link{ripple}}.}
 \item{lod.cutoff}{Only marker orders with LOD score (relative to the
   initial order) above this cutoff will be displayed.  For output of
   \code{\link{ripple}} in the case of minimization of the number
   of obligate crossovers, we double this argument and treat it as a
   cutoff for the number of obligate crossovers.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.ripple}, whose rows correspond to
  marker orders with likelihood (or number of obligate crossovers)
  within some cutoff of the initial order. If no marker order, other
  than the initial one, has likelihood within the specified range, the
  initial and next-best orders are returned.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{ % I don't run the following, since they are run in ripple.Rd
\dontrun{data(badorder)
rip1 <- ripple(badorder, 1, 7)
summary(rip1)

rip2 <- ripple(badorder, 1, 2, method="likelihood")
summary(rip2)

badorder <- switch.order(badorder, 1, rip2[2,])
}}

\seealso{ \code{\link{ripple}}, \code{\link{est.map}},
  \code{\link{est.rf}} }

\keyword{print}
