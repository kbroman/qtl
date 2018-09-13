\name{max.scanPhyloQTL}
\alias{max.scanPhyloQTL}

\title{Maximum peak in genome scan to map a QTL to a phylogenetic tree}

\description{
  Print the chromosome with the maximum LOD score across partitions, 
  from the results of \code{\link{scanPhyloQTL}}.
}

\usage{
\method{max}{scanPhyloQTL}(object, chr, format=c("postprob", "lod"), 
        \dots)
}

\arguments{
 \item{object}{An object output by the function
   \code{\link{scanPhyloQTL}}.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{format}{Indicates whether to provide LOD scores or approximate
   posterior probabilities; see the help file for \code{\link{summary.scanPhyloQTL}}.}
 \item{\dots}{Ignored at this point.}
}

\details{
  The output, and the use of the argument \code{format}, is as in
  \code{\link{summary.scanPhyloQTL}}. 
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
  \code{\link{summary.scanPhyloQTL}}, \code{\link{max.scanone}},
  \code{\link{inferredpartitions}}, 
  \code{\link{simPhyloQTL}} 
}

\references{
  Broman, K. W., Kim, S., An\'e, C. and Payseur, B. A.  Mapping
  quantitative trait loci to a phylogenetic tree.  In preparation.
}


\keyword{print}
