\name{markernames}
\alias{markernames}

\title{Pull out the marker names from a cross}

\description{
  Pull out the marker names from a cross object as one big vector.
}

\usage{markernames(cross, chr)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
}

\value{
  A vector of character strings (the marker names).
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
markernames(listeria, chr=5)
}

\seealso{ \code{\link{pull.map}}, \code{\link{phenames}}, \code{\link{chrnames}} }


\keyword{utilities}
