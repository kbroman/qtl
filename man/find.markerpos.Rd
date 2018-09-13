\name{find.markerpos}
\alias{find.markerpos}

\title{Find position of a marker}

\description{
  Find the chromosome and cM position of a set of genetic markers. 
}

\usage{find.markerpos(cross, marker)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{marker}{A vector of marker names.}
}

\value{
  A data frame with two columns: the chromosome and position of the
  markers.
}


\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
find.markerpos(hyper, "D4Mit164")
find.markerpos(hyper, c("D4Mit164", "D1Mit94"))
}

\seealso{ \code{\link{find.flanking}}, \code{\link{find.marker}},
\code{\link{find.pseudomarker}} }

\keyword{utilities}
