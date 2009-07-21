\name{subset.map}
\alias{subset.map}
\alias{[.map}

\title{Subsetting chromosomes for a genetic map}

\description{
  Pull out a specified set of chromosomes from a 
  \code{map} object.
}

\usage{
\method{subset}{map}(x, \dots)
\method{[}{map}(x, \dots)
}

\arguments{
 \item{x}{A list whose components are vectors of marker locations.}
 \item{\dots}{Vector of chromosome indices.}
}

\value{The input \code{map} object, but with only the specified subset
  of chromosomes.}


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(map10)
map10 <- subset(map10, chr=1:5)

# you can also use brackets
map10 <- map10[2:3]
}

\seealso{ \code{\link[qtl]{subset.cross}} }

\keyword{manip}
