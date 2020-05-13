\name{find_large_intervals}
\alias{find_large_intervals}

\title{Find large intervals in a map}

\description{
  Find large inter-marker intervals in a map.
}

\usage{
find_large_intervals(map, min_length=35)
}

\arguments{
 \item{map}{A list of numeric vectors; each component is a chromosome
   with the positions of markers on that chromosome. Can also
   be an object of class \code{cross}, in which case
 \code{\link{pull.map}} is used.}
 \item{min_length}{Minimum length of interval to be flagged.}
}


\value{
  Data frame with chromosome, left and right markers and interval length.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)

find_large_intervals(fake.f2, 30)
}

\seealso{ \code{\link{summaryMap}} }

\keyword{utilities}
