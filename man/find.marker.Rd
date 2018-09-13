\name{find.marker}
\alias{find.marker}

\title{Find marker closest to a specified position}

\description{
  Find the genetic marker closest to a specified position on a
  chromosome.
}

\usage{find.marker(cross, chr, pos, index)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{A vector of chromosome identifiers, or a single such.}
 \item{pos}{A vector of cM positions.}
 \item{index}{A vector of numeric indices of the markers within
   chromosomes.}
}

\value{
  A vector of marker names (of the same length as the input \code{pos}),
  corresponding to the markers nearest to the specified
  chromosomes/positions (if \code{pos} is specified) or to the input
  numeric indices (in \code{index} is specified).
}

\details{
  Provide one of \code{pos} or \code{index}.

  If the input \code{chr} has length one, it is expanded to the same
  length as the input \code{pos} or \code{index}.

  If \code{pos} is specified and multiple markers are exactly the same
  distance from the specified position, one is chosen at random from
  among those with the most genotype data.

  For a cross with sex-specific maps, positions specified by \code{pos}
  are assumed to correspond to the female genetic map.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
find.marker(listeria, 5, 28)
find.marker(listeria, 5, index=6)
find.marker(listeria, c(1, 5, 13), c(81, 28, 26))
}

\seealso{ \code{\link{find.flanking}}, \code{\link{plotPXG}},
  \code{\link{find.pseudomarker}}, \code{\link{effectplot}},
\code{\link{find.markerpos}} }

\keyword{utilities}
