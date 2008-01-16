\name{find.pseudomarker}
\alias{find.pseudomarker}

\title{Find the pseudomarker closest to a specified position}

\description{
  Find the pseudomarker closest to a specified position on a
  chromosome. 
}

\usage{find.pseudomarker(cross, chr, pos, where=c("draws", "prob"))}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector of chromosome identifiers, or a single such.}
 \item{pos}{A vector of cM positions.}
 \item{where}{Indicates whether to look in the \code{draws} or
 \code{prob} components of the input cross.}
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

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(listeria)
listeria <- calc.genoprob(listeria, step=2.5)
find.pseudomarker(listeria, 5, 28, "prob")
find.pseudomarker(listeria, c(1, 5, 13), c(81, 28, 26), "prob")
}

\seealso{ \code{\link[qtl]{find.flanking}}, \code{\link[qtl]{plot.pxg}},
\code{\link[qtl]{effectplot}}, \code{\link[qtl]{find.marker}},
\code{\link[qtl]{find.markerpos}}  }

\keyword{utilities}
