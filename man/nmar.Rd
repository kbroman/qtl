\name{nmar}
\alias{nmar}

\title{Determine the numbers of markers on each chromosome}

\description{
  Determine the number of markers on each chromosome in a cross
  or map object.
}

\usage{
nmar(object)
}

\arguments{
 \item{object}{An object of class \code{cross} (see
   \code{\link{read.cross}} for details) or \code{map} (see
 \code{\link{sim.map}} for details).}
}

\value{
  A vector with the numbers of markers on each chromosome in the input.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
nmar(fake.f2)
map <- pull.map(fake.f2)
nmar(map)
}

\seealso{ \code{\link{read.cross}}, \code{\link{plot.cross}},
  \code{\link{summary.cross}},
  \code{\link{nind}},
  \code{\link{nchr}},
  \code{\link{totmar}},
  \code{\link{nphe}} }

\keyword{print}
