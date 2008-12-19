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
   \code{\link[qtl]{read.cross}} for details) or \code{map} (see
 \code{\link[qtl]{sim.map}} for details).}
}

\value{
  The number of chromosomes in the input.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
nchr(fake.f2)
map <- pull.map(fake.f2)
nchr(map)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}},
  \code{\link[qtl]{summary.cross}},
  \code{\link[qtl]{nind}},
  \code{\link[qtl]{totmar}},
  \code{\link[qtl]{nmar}},
  \code{\link[qtl]{nphe}} }

\keyword{print}
