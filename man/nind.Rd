\name{nind}
\alias{nind}

\title{Determine the number of individuals QTL experiment}

\description{
  Determine the number of individuals in cross object.
}

\usage{
nind(object)
}

\arguments{
 \item{object}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  The number of individuals in the input cross object.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
nind(fake.f2)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}},
  \code{\link[qtl]{summary.cross}},
  \code{\link[qtl]{nmar}},
  \code{\link[qtl]{nchr}},
  \code{\link[qtl]{totmar}},
  \code{\link[qtl]{nphe}} }

\keyword{print}
