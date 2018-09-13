\name{nchr}
\alias{nchr}

\title{Determine the number of chromosomes}

\description{
  Determine the number of chromosomes in a cross
  or map object.
}

\usage{
nchr(object)
}

\arguments{
 \item{object}{An object of class \code{cross} (see
   \code{\link{read.cross}} for details) or \code{map} (see
 \code{\link{sim.map}} for details).}
}

\value{
  The number of chromosomes in the input.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
nchr(fake.f2)
map <- pull.map(fake.f2)
nchr(map)
}

\seealso{ \code{\link{read.cross}}, \code{\link{plot.cross}},
  \code{\link{summary.cross}},
  \code{\link{nind}},
  \code{\link{totmar}},
  \code{\link{nmar}},
  \code{\link{nphe}} }

\keyword{print}
