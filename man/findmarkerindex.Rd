\name{find.markerindex}
\alias{find.markerindex}

\title{Determine the numeric index for a marker}

\description{
  Determine the numeric index for a marker in a cross object, when all
  markers on all chromosomes are pasted together.
}

\usage{find.markerindex(cross, name)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{name}{A vector of marker names.}
}

\value{
  A vector of numeric indices, from 1, 2, \dots, \code{totmar(cross)},
  with \code{NA} for markers not found. 
}


\author{Danny Arends; Karl W Broman \email{broman@wisc.edu} }

\examples{
data(hyper)
mar <- find.marker(hyper, 4, 30)
find.markerindex(hyper, mar)
}

\seealso{ \code{\link{find.markerpos}} }

\keyword{utilities}
