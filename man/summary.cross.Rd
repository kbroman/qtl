\name{summary.cross}
\alias{summary.cross}

\title{Print summary of QTL experiment}

\description{
  Print summary information about a \code{cross} object.
}

\usage{
\method{summary}{cross}(object, \dots)
}

\arguments{
 \item{object}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.cross}
  containing a variety of summary information about the cross (this is
  generally printed automatically).
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
summary(fake.f2)
}

\seealso{ \code{\link{read.cross}}, \code{\link{plot.cross}},
  \code{\link{nind}},
  \code{\link{nmar}},
  \code{\link{nchr}},
  \code{\link{totmar}},
  \code{\link{nphe}} }

\keyword{print}
