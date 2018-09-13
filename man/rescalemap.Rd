\name{rescalemap}
\alias{rescalemap}

\title{Rescale genetic maps}

\description{
  Rescale a genetic map by multiplying all positions by a constant
}

\usage{
rescalemap(object, scale=1e-6)
}
\arguments{
 \item{object}{An object of class \code{cross} (see
   \code{\link{read.cross}} for details) or \code{map} (see
   \code{\link{sim.map}} for details).}
 \item{scale}{Scale factor by which all positions will be multiplied.}
}

\details{
  This function is included particularly for the case that map positions
  in a cross object were provided in basepairs and one wishes to quickly
  convert them to Mbp or some other approximation of cM distances.  (In
  the mouse, 1 cM is approximation 2 Mbp, so one might use
  \code{scale=5e-7} in this function.)
}

\value{
If the input is a \code{map} object, a \code{map} object is returned; if
the input is a \code{cross} object, a \code{cross} object is returned.
In either case, the positions of markers are simply multiplied by
\code{scale}.  

}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
rescaled <- rescalemap(hyper, scale=2)
plotMap(hyper, rescaled)
}


\seealso{ \code{\link{replace.map}},
  \code{\link{est.map}} }

\keyword{utilities}
