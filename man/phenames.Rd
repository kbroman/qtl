\name{phenames}
\alias{phenames}

\title{Pull out the phenotypes names from a cross}

\description{
  Pull out the phenotype names from a cross object as a vector.
}

\usage{phenames(cross)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
}

\value{
  A vector of character strings (the phenotype names).
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
phenames(listeria)
}

\seealso{ \code{\link{markernames}}, \code{\link{chrnames}} }


\keyword{utilities}
