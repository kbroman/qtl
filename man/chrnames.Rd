\name{chrnames}
\alias{chrnames}

\title{Pull out the chromosome names from a cross}

\description{
  Pull out the chromosome names from a cross object as one big vector.
}

\usage{chrnames(cross)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
}

\value{
  A vector of character strings (the chromosome names).
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
chrnames(listeria)
}

\seealso{ \code{\link{markernames}}, \code{\link{phenames}} }


\keyword{utilities}
