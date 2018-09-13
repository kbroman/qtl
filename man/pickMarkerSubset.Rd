\name{pickMarkerSubset}
\alias{pickMarkerSubset}

\title{Identify the largest subset of markers that are some distance apart}

\description{
  Identify the largest subset of markers for which no two adjacent
  markers are separated by less than some specified distance; if weights
  are provided, find the marker subset for which the sum of the weights
  is maximized.
}

\usage{
pickMarkerSubset(locations, min.distance, weights)
}

\arguments{
 \item{locations}{A vector of marker locations.}
 \item{min.distance}{Minimum distance between adjacent markers in the
 chosen subset.}
 \item{weights}{(Optional) vector of weights for the markers.  If
 missing, we take \code{weights == 1}.}
}

\value{A vector of marker names.}


\details{
Let \eqn{d_i}{d[i]} be
the location of marker \eqn{i}, for \eqn{i \in 1, \dots, M}{i in 1,
\dots, M}.  We use the dynamic programming algorithm of Broman and Weber
(1999) to identify the subset of markers \eqn{i_1, \dots, i_k}{i[1], \dots,
i[k]} for which \eqn{d_{i_{j+1}} - d_{i_j} \le}{d(i[j+1]) - d(i[j]) <=}
\code{min.distance} and \eqn{\sum w_{i_j}}{sum w(i[j])} is maximized.

If there are multiple optimal subsets, we pick one at random.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\references{
  Broman,  K. W. and Weber,  J. L. (1999) Method for constructing
  confidently ordered linkage maps. \emph{Genet. Epidemiol.}, \bold{16},
  337--343.
}

\examples{
data(hyper)

# subset of markers on chr 4 spaced >= 5 cM
pickMarkerSubset(pull.map(hyper)[[4]], 5)

# no. missing genotypes at each chr 4 marker
n.missing <- nmissing(subset(hyper, chr=4), what="mar")

# weight by -log(prop'n missing), but don't let 0 missing go to +Inf
wts <- -log( (n.missing+1) / (nind(hyper)+1) )

# subset of markers on chr 4 spaced >= 5 cM, with weights = -log(prop'n missing)
pickMarkerSubset(pull.map(hyper)[[4]], 5, wts)
}

\seealso{ \code{\link{drop.markers}}, \code{\link{pull.markers}} }

\keyword{manip}
