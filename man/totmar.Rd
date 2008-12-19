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
   \code{\link[qtl]{read.cross}} for details) or \code{map} (see
 \code{\link[qtl]{sim.map}} for details).}
}

\value{
  The total number of markers in the input.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
totmar(fake.f2)
map <- pull.map(fake.f2)
totmar(map)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}},
  \code{\link[qtl]{summary.cross}},
  \code{\link[qtl]{nind}},
  \code{\link[qtl]{nchr}},
  \code{\link[qtl]{nmar}},
  \code{\link[qtl]{nphe}} }

\keyword{print}
