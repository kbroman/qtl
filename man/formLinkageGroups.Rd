\name{formLinkageGroups}
\alias{formLinkageGroups}

\title{Partition markers into linkage groups}

\description{
  Use pairwise linkage information between markers (as calculated by
  \code{\link{est.rf}} to partition markers into linkage groups.
}

\usage{
formLinkageGroups(cross, max.rf=0.25, min.lod=3, reorgMarkers=FALSE,
                  verbose=FALSE)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{max.rf}{Maximum recombination fraction for placing two markers in
   the same linkage group (see Details).}
 \item{min.lod}{Minimum LOD score for placing two markers in the same
 linkage group (see Details).}
 \item{reorgMarkers}{If TRUE, the output is a cross object, like the
 input, but with the markers organized into the inferred linkage
 groups.  If FALSE, the output is a table indicating the initial
 chromosome assignments and the inferred linkage group partitions.}
 \item{verbose}{If TRUE, display information about the progress of the
 calculations.}
}

\value{
  If \code{reorgMarkers=FALSE} (the default), the output is a data frame
  with rows corresponding to the markers and with two columns: the
  initial chromosome assignment and the inferred linkage group.  Linkage
  groups are ordered by the number of markers they contain (from largest
  to smallest).

  If \code{reorgMarkers=TRUE}, the output is a cross object, like the
  input, but with the markers reorganized into the inferred linkage
  groups.  The marker order and marker positions within the linkage
  groups are arbitrary.
}

\details{
  Two markers are placed in the same linkage group if the estimated
  recombination fraction between them is \eqn{\le}{<=} \code{max.rf} and
  the LOD score (for the test of the rec. frac. = 1/2) is \eqn{\ge}{>=}
  \code{min.lod}.  The transitive property (if A is linked to B and B is
  linked to C then A is linked to C) is used to close the groups.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
listeria <- est.rf(listeria)
result <- formLinkageGroups(listeria)
tab <- table(result[,1], result[,2])
apply(tab, 1, function(a) sum(a!=0))
apply(tab, 2, function(a) sum(a!=0))
}

\seealso{ \code{\link{est.rf}}, \code{\link{orderMarkers}} }


\keyword{utilities}
