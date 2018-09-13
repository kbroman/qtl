\name{find.pseudomarker}
\alias{find.pseudomarker}

\title{Find the pseudomarker closest to a specified position}

\description{
  Find the pseudomarker closest to a specified position on a
  chromosome. 
}

\usage{find.pseudomarker(cross, chr, pos, where=c("draws", "prob"), addchr=TRUE)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{A vector of chromosome identifiers, or a single such.}
 \item{pos}{A vector of cM positions.}
 \item{where}{Indicates whether to look in the \code{draws} or
 \code{prob} components of the input cross.}
 \item{addchr}{If TRUE, include something like \code{"c5."} at the
 beginning of the names of non-pseudomarker locations, as in the output
 of \code{\link{scanone}}; if FALSE, don't include this sort of string,
 as in the genotype probabilities from \code{\link{calc.genoprob}}.}
}

\value{
  A vector of pseudomarker names (of the same length as the input \code{pos}),
  corresponding to the markers nearest to the specified chromosomes/positions. 
}

\details{
  If the input \code{chr} has length one, it is expanded to the same
  length as the input \code{pos}.

  If multiple markers are exactly the same distance from the specified
  position, one is chosen at random from among those with the most
  genotype data.

  For a cross with sex-specific maps, the input positions are assumed to
  correspond to the female genetic map.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
listeria <- calc.genoprob(listeria, step=2.5)
find.pseudomarker(listeria, 5, 28, "prob")
find.pseudomarker(listeria, c(1, 5, 13), c(81, 28, 26), "prob")
}

\seealso{ \code{\link{find.flanking}}, \code{\link{plotPXG}},
\code{\link{effectplot}}, \code{\link{find.marker}},
\code{\link{find.markerpos}}  }

\keyword{utilities}
