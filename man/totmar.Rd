\name{totmar}
\alias{totmar}

\title{Determine the total number of markers}

\description{
  Determine the total number of markers in a cross
  or map object.
}

\usage{
totmar(object)
}

\arguments{
 \item{object}{An object of class \code{cross} (see
   \code{\link{read.cross}} for details) or \code{map} (see
 \code{\link{sim.map}} for details).}
}

\value{
  The total number of markers in the input.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
totmar(fake.f2)
map <- pull.map(fake.f2)
totmar(map)
}

\seealso{ \code{\link{read.cross}}, \code{\link{plot.cross}},
  \code{\link{summary.cross}},
  \code{\link{nind}},
  \code{\link{nchr}},
  \code{\link{nmar}},
  \code{\link{nphe}} }

\keyword{print}
