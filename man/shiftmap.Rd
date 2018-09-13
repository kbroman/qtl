\name{shiftmap}
\alias{shiftmap}

\title{Shift starting points in genetic maps}

\description{
  Shift starting points in a genetic map to a set of defined positions
}

\usage{
shiftmap(object, offset=0)
}
\arguments{
 \item{object}{An object of class \code{cross} (see
   \code{\link{read.cross}} for details) or \code{map} (see
   \code{\link{sim.map}} for details).}
 \item{offset}{Defines the starting position
  for each chromosome.  This should be a
  single value (to be used for all chromosomes) or a vector with length
  equal to the number of chromosomes, defining individual starting
  positions for each chromosome.  For a sex-specific map (as in a 4-way
  cross), we use the same offset for both the male and female maps.}
}


\value{
If the input is a \code{map} object, a \code{map} object is returned; if
the input is a \code{cross} object, a \code{cross} object is returned.
In either case, the positions of markers are shifted so that the
starting positions are as in \code{offset}.  

}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
shiftedhyper <- shiftmap(hyper, offset=0)
par(mfrow=c(1,2))
plotMap(hyper, shift=FALSE, alternate.chrid=TRUE)
plotMap(shiftedhyper, shift=FALSE, alternate.chrid=TRUE)
}


\seealso{ \code{\link{replace.map}},
  \code{\link{est.map}} }

\keyword{utilities}
