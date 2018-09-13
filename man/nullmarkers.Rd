\name{nullmarkers}
\alias{nullmarkers}

\title{Identify markers without any genotype data}

\description{
  Identify markers in a cross that have no genotype data.
}

\usage{
nullmarkers(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
}

\value{
  Marker names (a vector of character strings) with no genotype data.
}

\examples{
# one marker with no data
data(hyper)
nullmarkers(hyper)

# nothing in listeria
data(listeria)
nullmarkers(listeria)
}

\seealso{ \code{\link{drop.nullmarkers}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{utilities}
