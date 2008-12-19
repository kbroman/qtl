\name{nphe}
\alias{nphe}

\title{Determine the number of phenotypes QTL experiment}

\description{
  Determine the number of phenotypes in cross object.
}

\usage{
nphe(object)
}

\arguments{
 \item{object}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  The number of phenotypes in the input cross object.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
nphe(fake.f2)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}},
  \code{\link[qtl]{summary.cross}},
  \code{\link[qtl]{nmar}},
  \code{\link[qtl]{nchr}},
  \code{\link[qtl]{totmar}},
  \code{\link[qtl]{nind}} }

\keyword{print}
