\name{markerlrt}
\alias{markerlrt}

\title{General likelihood ratio test for association between marker pairs}

\description{
  Calculate a LOD score for a general likelihood ratio test for each
  pair of markers, to assess their association.
}

\usage{
markerlrt(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
}

\value{
  The input \code{cross} object is returned with a component, \code{rf},
  added.  This is a matrix of size (tot.mar x tot.mar).  The diagonal
  contains the number of typed meioses per marker, the upper and lower triangles
  each contain the LOD scores.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(badorder)
badorder <- markerlrt(badorder)
plotRF(badorder)
}

\seealso{ \code{\link{plotRF}}, \code{\link{est.rf}},
  \code{\link{badorder}} }

\keyword{utilities}
