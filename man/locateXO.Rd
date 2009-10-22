\name{locateXO}
\alias{locateXO}

\title{Estimate locations of crossovers}

\description{
  Estimate the locations of crossovers for each individual on a given
  chromosome.
}

\usage{
locateXO(cross, chr)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Chromosome to investigate (if unspecified, the first
   chromosome is considered).
   This should be a character string referring to a chromosome
   by name; numeric values are converted to strings.}
}

\details{
  For each individual we detemine the locations of obligate crossovers,
  and estimate their location to be at the midpoint between the nearest
  flanking typed markers.

  The function currently only works for a backcross, intercross, or
  recombinant inbred line.
}

\value{
  A list with one component per individual.  Each component is either
  NULL or is a numeric vector with the estimated crossover locations.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(hyper)
xoloc <- locateXO(hyper, chr=4)
table(sapply(xoloc, length))
}

\seealso{ \code{\link{countXO}} }

\keyword{utilities}
