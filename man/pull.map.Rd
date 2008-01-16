\name{pull.map}
\alias{pull.map}

\title{Pull out the genetic map from a cross}

\description{
  Pull out the map portion of a cross object.
}

\usage{pull.map(cross, chr)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{An optional vector specifying which chromosomes to keep or
   discard. This may be a logical vector, a numeric vector, or a vector
   of character strings.} 
}

\value{
  The genetic map: a list with each component containing the marker
  positions (in cM) for a chromosome.  Each component has class
  \code{A} or \code{X} according to whether it is an autosome or the X
  chromosome.  The components are either vectors of marker positions or,
  for a sex-specific map, 2-row matrices containing the female and male
  marker locations.  The map itself is given class \code{map}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
map <- pull.map(fake.f2)
plot(map)
}

\seealso{ \code{\link[qtl]{replace.map}}, \code{\link[qtl]{plot.map}} }


\keyword{utilities}
