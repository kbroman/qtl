\name{orderMarkers}
\alias{orderMarkers}

\title{Find an initial order for markers within chromosomes}

\description{
  Establish initial orders for markers within chromosomes by a greedy
  algorithm, adding one marker at a time with locations of previous
  markers fixed, in the position giving the miniminum number of obligate
  crossovers.
}

\usage{
orderMarkers(cross, chr, window=7, use.ripple=TRUE, error.prob=0.0001,
             map.function=c("haldane","kosambi","c-f","morgan"),
             maxit=4000, tol=1e-4, sex.sp=TRUE, verbose=FALSE)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
   This should be a vector of character strings referring to chromosomes
   by name; numeric values are converted to strings.  Refer to
   chromosomes with a preceding \code{-} to have all chromosomes but
   those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{window}{If \code{use.ripple=TRUE}, this indicates the number of
   markers to include in the sliding window of permuted markers.  Larger
   numbers result in the comparison of a greater number of marker
   orders, but will require a considerable increase in computation
   time.}
 \item{use.ripple}{If TRUE, the initial order is refined by a call to
 the function \code{\link{ripple}}.}
 \item{error.prob}{Assumed genotyping error rate used in the final
   estimated map.}
 \item{map.function}{Indicates the map function to use in the final
   estimated map.}
 \item{maxit}{Maximum number of EM iterations to perform in the final
   estimated map.}
 \item{tol}{Tolerance for determining convergence in the final estimated
   map.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps in the
   final estimated map; this is used only for the 4-way cross.}
 \item{verbose}{If TRUE, information about the progress of the
   calculations is displayed; if > 1, even more information is given.}
}

\value{
  The output is a cross object, as in the input, with orders of markers
  on selected chromosomes revised.
}

\details{
  Markers within a linkage group are considered in order of decreasing
  number of genotyped individuals.  The first two markers are placed in
  an arbitrary order.  Additional markers are considered one at a time,
  and each possible placement of a marker is compared (with the order of
  the previously placed markers taken as fixed) via the number of
  obligate crossovers (that is, the minimal number of crossovers that
  would explain the observed data).  The marker is placed in the
  position giving the minimal number of obligate crossovers.  If
  multiple positions give the same number of obligate crossovers, a
  single location (among those positions) is chosen at random.

  If \code{use.ripple=TRUE}, the final order is passed to
  \code{\link{ripple}} with \code{method="countxo"} to refine the
  marker order.  If \code{use.ripple=TRUE} and the number of markers on
  a chromosome is \eqn{\le}{<=} the argument \code{window}, the initial
  greedy algorithm is skipped and all possible marker orders are
  compared via \code{\link{ripple}}.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
pull.map(listeria, chr=3)
revcross <- orderMarkers(listeria, chr=3, use.ripple=FALSE)
pull.map(revcross, chr=3)
}

\seealso{ \code{\link{formLinkageGroups}},
  \code{\link{ripple}}, \code{\link{est.map}}, \code{\link{countXO}} }

\keyword{utilities}
