\name{ripple}
\alias{ripple}

\title{Compare marker orders}

\description{
  Investigate different marker orders for a given chromosome, comparing
  all possible permutations of a sliding window of markers.
}

\usage{
ripple(cross, chr, window=4, method=c("countxo","likelihood"),
       error.prob=0.0001, map.function=c("haldane","kosambi","c-f","morgan"),
       maxit=4000, tol=1e-6, sex.sp=TRUE, verbose=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{The chromosome to investigate.  Only one chromosome is
  allowed. (This should be a character string referring to the
  chromosomes by name.)}
 \item{window}{Number of markers to include in the sliding window of
   permuted markers.  Larger numbers result in the comparison of a
   greater number of marker orders, but will require a considerable
   increase in computation time.}
 \item{method}{Indicates whether to compare orders by counting the
   number of obligate crossovers, or by a likelihood analysis.} 
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
 \item{maxit}{Maximum number of EM iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps; this is 
   used only for the 4-way cross.}
 \item{verbose}{If TRUE, information about the number of orders (and, if
   \code{method="likelihood"}, about progress) are printed.}
}

\details{
  For \code{method="likelihood"}, calculations are done by first
  constructing a matrix of marker orders and then making repeated calls
  to the R function \code{\link[qtl]{est.map}}.  Of course, it would be
  faster to do everything within C, but this was a lot easier to code.

  For \code{method="countxo"}, calculations are done within C.  
}

\value{
  A matrix, given class \code{"ripple"}; the first set of columns are
  marker indices describing the order.  In the case of
  \code{method="countxo"}, the last column is the number of obligate
  crossovers for each particular order.  In the case of
  \code{method="likelihood"}, the last two columns are LOD scores (log
  base 10 likelihood ratios) comparing each order to the initial order
  and the estimated chromosome length for the given order.  Positive LOD
  scores indicate that the alternate order has more support than the
  original.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(badorder)
rip1 <- ripple(badorder, chr=1, window=3)
summary(rip1)

\dontrun{
rip2 <- ripple(badorder, chr=1, window=2, method="likelihood")
summary(rip2)
}

badorder <- switch.order(badorder, 1, rip1[2,])
}

\seealso{ \code{\link[qtl]{summary.ripple}}, \code{\link[qtl]{switch.order}},
  \code{\link[qtl]{est.map}}, \code{\link[qtl]{est.rf}} }

\keyword{utilities}
