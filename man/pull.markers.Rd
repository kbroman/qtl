\name{pull.markers}
\alias{pull.markers}

\title{Drop all but a selected set of markers}

\description{
  Drop all but a selected set of markers from the data matrices and genetic maps.
}

\usage{
pull.markers(cross, markers)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{markers}{A character vector of marker names.}
}

\value{
  The input object, with any markers not specified in the vector \code{markers} removed
  from the genotype data matrices, genetic maps, and, if applicable, any
  derived data (such as produced by \code{\link{calc.genoprob}}).
  (It might be a good idea to re-derive such things after using this
  function.)
}

\examples{
data(listeria)
listeria2 <- pull.markers(listeria, c("D10M44","D1M3","D1M75"))
}

\seealso{ \code{\link{drop.nullmarkers}}, \code{\link{drop.markers}},
  \code{\link{geno.table}},
  \code{\link{clean.cross}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{manip}
