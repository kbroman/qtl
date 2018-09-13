\name{chrlen}
\alias{chrlen}

\title{Chromosome lengths in QTL experiment}

\description{
  Obtain the chromosome lengths in a \code{cross} or \code{map} object.
}

\usage{
chrlen(object)
}

\arguments{
 \item{object}{An object of class \code{map} or of class \code{cross}.}
}

\value{
  Returns a vector of chromosome lengths.  If the cross has sex-specific
  maps, it returns a 2-row matrix with the two lengths for each chromosome.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
chrlen(fake.f2)

map <- pull.map(fake.f2)
chrlen(map)
}

\seealso{ \code{\link{summaryMap}}, \code{\link{pull.map}},
  \code{\link{summary.cross}} } 

\keyword{print}
