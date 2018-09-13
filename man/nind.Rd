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
   \code{\link{read.cross}} for details.}
}

\value{
  The number of individuals in the input cross object.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
nind(fake.f2)
}

\seealso{ \code{\link{read.cross}}, \code{\link{plot.cross}},
  \code{\link{summary.cross}},
  \code{\link{nmar}},
  \code{\link{nchr}},
  \code{\link{totmar}},
  \code{\link{nphe}} }

\keyword{print}
