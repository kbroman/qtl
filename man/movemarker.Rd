\name{movemarker}
\alias{movemarker}
\title{Move a marker to a new chromosome}

\description{
  Move a specified marker to a different chromosome.
}

\usage{movemarker(cross, marker, newchr, newpos)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{marker}{The name of the marker to be moved (a character string).}
 \item{newchr}{The chromosome to which the marker should be moved.}
 \item{newpos}{The position (in cM) at which the marker should be
   placed.  If missing, the marker is placed at the end of the
   chromosome.}
}

\value{
  The input \code{cross} object, but with the specified marker moved to
  the specified chromosome.

  All intermediate calculations (such as the results of
  \code{\link{calc.genoprob}} and \code{\link{est.rf}}) are
  removed.
}


\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(badorder)
badorder <- movemarker(badorder, "D2M937", 3, 48.15)
badorder <- movemarker(badorder, "D3M160", 2, 28.83)
}

\seealso{ \code{\link{switch.order}} }

\keyword{manip}
