\name{convert.map}
\alias{convert.map}

\title{Change map function for a genetic map}

\description{
  Convert a genetic map from using one map function to another.
}

\usage{
\method{convert}{map}(object, old.map.function=c("haldane", "kosambi", "c-f", "morgan"),
         new.map.function=c("haldane", "kosambi", "c-f", "morgan"), \dots)
}

\arguments{
 \item{object}{A genetic map object, of class \code{"map"}: A list whose
 components are vectors of marker locations.} 
 \item{old.map.function}{The map function used in forming the map in
   \code{object}.}
 \item{new.map.function}{The new map function to be used.}
 \item{\dots}{Ignored at this point.}
}

\value{
  The same as the input, but with inter-marker distances changed to
  reflect a different map function.
}

\details{
  The location of the first marker on each chromosome is left
  unchanged.  Inter-marker distances are converted to recombination
  fractions with the inverse of the \code{old.map.function}, and then
  back to distances with the \code{new.map.function}.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
map <- pull.map(listeria)
map <- convert(map, "haldane", "kosambi")
listeria <- replace.map(listeria, map)
}

\seealso{ \code{\link{est.map}}, \code{\link{replace.map}} }

\keyword{manip}
