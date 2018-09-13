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
   \code{\link{read.cross}} for details.}
}

\value{
  The number of phenotypes in the input cross object.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
nphe(fake.f2)
}

\seealso{ \code{\link{read.cross}}, \code{\link{plot.cross}},
  \code{\link{summary.cross}},
  \code{\link{nmar}},
  \code{\link{nchr}},
  \code{\link{totmar}},
  \code{\link{nind}} }

\keyword{print}
