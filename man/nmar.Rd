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
   \code{\link[qtl]{read.cross}} for details) or \code{map} (see
 \code{\link[qtl]{sim.map}} for details).}
}

\value{
  A vector with the numbers of markers on each chromosome in the input.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
nmar(fake.f2)
map <- pull.map(fake.f2)
nmar(map)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}},
  \code{\link[qtl]{summary.cross}},
  \code{\link[qtl]{nind}},
  \code{\link[qtl]{nchr}},
  \code{\link[qtl]{totmar}},
  \code{\link[qtl]{nphe}} }

\keyword{print}
