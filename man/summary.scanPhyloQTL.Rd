\name{summary.scanPhyloQTL}
\alias{summary.scanPhyloQTL}

\title{Summarize the results a genome scan to map a QTL to a phylogenetic tree}

\description{
  Print the maximum LOD scores for each partition on each chromosome,
  from the results of \code{\link{scanPhyloQTL}}.
}

\usage{
\method{summary}{scanPhyloQTL}(object, format=c("postprob", "lod"),
        threshold, \dots)
}

\arguments{
 \item{object}{An object output by the function
   \code{\link{scanPhyloQTL}}.}
 \item{format}{Indicates whether to provide LOD scores or approximate
   posterior probabilities; see Details below.}
 \item{threshold}{A threshold determining which chromosomes should be
   output; see Details below.}
 \item{\dots}{Ignored at this point.}
}

\details{
  This function is used to report chromosomes deemed interesting from a one-QTL
  genome scan to map QTL to a phylogenetic tree (by \code{\link{scanPhyloQTL}}).

  For \code{format="lod"}, the output contains the maximum LOD score for
  each partition on each chromosome (which do not necessarily occur at
  the same position).  The position corresponds to the peak location for
  the partition with the largest LOD score on that chromosome.  The
  last column is the overall maximum LOD (across partitions) on that
  chromosome.  The second-to-last column is the inferred partition
  (i.e., that with the largest LOD
  score.  The third-to-last column is the difference between the LOD score for
  the best partition and that for the second-best.

  For \code{format="postprob"}, the final column contains the maximum
  LOD score across partitions.  But instead of providing the LOD
  scores for each partition, these are converted to approximate
  posterior probabilities under the assumption of a single diallelic QTL
  on that chromosome:  on each chromosome, we take
  \eqn{10^{LOD}}{10^LOD} for the partitions and rescale them to sum to 1.

  The \code{threshold} argument is applied to the last column (the
  maximum LOD score across partitions).
}

\value{
  An object of class \code{summary.scanPhyloQTL}, to be printed by
  \code{print.summary.scanPhyloQTL}.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
\dontrun{
# example map; drop X chromosome
data(map10)
map10 <- map10[1:19]

# simulate data
x <- simPhyloQTL(4, partition="AB|CD", crosses=c("AB", "AC", "AD"),
                 map=map10, n.ind=150,
                 model=c(1, 50, 0.5, 0))

# run calc.genoprob on each cross
x <- lapply(x, calc.genoprob, step=2)

# scan genome, at each position trying all possible partitions
out <- scanPhyloQTL(x, method="hk")

# maximum peak
max(out, format="lod")

# approximate posterior probabilities at peak
max(out, format="postprob")

# all peaks above a threshold for LOD(best) - LOD(2nd best)
summary(out, threshold=1, format="lod")

# all peaks above a threshold for LOD(best), showing approx post'r prob
summary(out, format="postprob", threshold=3)

# plot of results
plot(out)
}
}

\seealso{ \code{\link{scanPhyloQTL}}, \code{\link{plot.scanPhyloQTL}},
  \code{\link{max.scanPhyloQTL}}, \code{\link{summary.scanone}},
  \code{\link{inferredpartitions}},
  \code{\link{simPhyloQTL}}
}

\references{
  Broman, K. W., Kim, S., An\'e, C. and Payseur, B. A.  Mapping
  quantitative trait loci to a phylogenetic tree.  In preparation.
}


\keyword{print}
