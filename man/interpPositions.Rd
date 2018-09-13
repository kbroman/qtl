\name{interpPositions}
\alias{interpPositions}

\title{Interpolate positions from one map to another}

\description{
  On the basis of a pair of marker maps with common markers, take
  positions along one map and interpolate (or, past the terminal markers
  on a chromosome, extrapolate) their positions on the second map.
}

\usage{
interpPositions(oldpositions, oldmap, newmap)
}

\arguments{
 \item{oldpositions}{A data frame with two columns: \code{chr}
 (chromosome identifiers) and
 \code{pos} (positions, along \code{oldmap}).}
 \item{oldmap}{An object of class \code{"map"}; see
 \code{\link{sim.map}} for details.}
 \item{newmap}{An object of class \code{"map"}, with the same
 chromosomes and markers as \code{oldmap}.}
}

\value{The input data frame, \code{oldpositions}, with an additional
  column \code{newpos} with the interpolated positions along
  \code{newmap}.} 


\details{
In this explanation, take \code{oldmap} and \code{newmap} to be the
physical and genetic maps, respectively.  

We use linear interpolation within each interval, assuming a constant
recombination rate within the interval.  Past the terminal
markers, we use linear extrapolation, using the chromosome-wide average
recombination rate.
}

\author{Karl W Broman, \email{broman@wisc.edu} }


\examples{
data(hyper)

# hyper genetic map
gmap <- pull.map(hyper)

# a fake physical map, with each chromosome starting at 0.
pmap <- shiftmap(rescalemap(gmap, 2))

# positions on pmap to determine location on gmap
tofind <- data.frame(chr=c(1, 5, 17, "X"), pos=c(220, 20, 105, 10))
rownames(tofind) <- paste("loc", 1:nrow(tofind), sep="")

interpPositions(tofind, pmap, gmap)
}

\seealso{ \code{\link{shiftmap}}, \code{\link{rescalemap}}, \code{\link{pull.map}} }

\keyword{manip}
