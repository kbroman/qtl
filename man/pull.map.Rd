\name{pull.map}
\alias{pull.map}

\title{Pull out the genetic map from a cross}

\description{
  Pull out the map portion of a cross object.
}

\usage{pull.map(cross, chr, as.table=FALSE)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{as.table}{If TRUE, return the genetic map as a table with
 chromosome assignments and marker names.  If FALSE, return the map as a
 \code{"map"} object.}
}

\value{
  The genetic map: a list with each component containing the marker
  positions (in cM) for a chromosome.  Each component has class
  \code{A} or \code{X} according to whether it is an autosome or the X
  chromosome.  The components are either vectors of marker positions or,
  for a sex-specific map, 2-row matrices containing the female and male
  marker locations.  The map itself is given class \code{map}.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
map <- pull.map(fake.f2)
plot(map)
}

\seealso{ \code{\link{replace.map}}, \code{\link{plotMap}}, \code{\link{map2table}} }


\keyword{utilities}
