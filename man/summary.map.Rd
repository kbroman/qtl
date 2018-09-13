\name{summaryMap}
\alias{summary.map}
\alias{summaryMap}

\title{Print summary of a genetic map}

\description{
  Print summary information about a \code{map} object.
}

\usage{
\method{summary}{map}(object, \dots)
summaryMap(object, \dots)
}

\arguments{
 \item{object}{An object of class \code{map}, which is a list of
   vectors (or, for a sex-specific map, 2-row matrices), each specifying
   the locations of the markers.  The object can also be of class
   \code{cross}, in which case the function \code{\link{pull.map}}
   is used to extract the genetic map from the object.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.map}, which is just a data.frame 
  containing the number of markers, length, the average inter-marker
  spacing, and the maximum distance between markers, for each chromosome
  and overall.  An attribute \code{sexsp} indicates whether the map was
  sex-specific.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(map10)
summary(map10)
}

\seealso{ \code{\link{chrlen}}, \code{\link{pull.map}},
  \code{\link{summary.cross}} } 

\keyword{print}
