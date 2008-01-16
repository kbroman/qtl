\name{find.flanking}
\alias{find.flanking}

\title{Find flanking markers for a specified position}

\description{
  Find the genetic markers flanking a specified position on a
  chromosome, as well as the marker that is closest to the specified
  position.
}

\usage{find.flanking(cross, chr, pos)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector of chromosome identifiers, or a single such.}
 \item{pos}{A vector of cM positions.}
}

\value{
  A data.frame, each row corresponding to one of the input positions.
  The first column contains the left-flanking markers, the second column
  contains the right-flanking markers, and the third column contains the
  markers closest to the specified positions.
}

%\details{
%}

\author{Brian Yandell}

\examples{
data(listeria)
find.flanking(listeria, 5, 28)
find.flanking(listeria, c(1, 5, 13), c(81, 28, 26))
}

\seealso{ \code{\link[qtl]{find.marker}}, \code{\link[qtl]{plot.pxg}},
\code{\link[qtl]{find.markerpos}}, \code{\link[qtl]{find.pseudomarker}} }

\keyword{utilities}
