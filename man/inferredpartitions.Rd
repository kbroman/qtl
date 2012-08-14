\name{inferredpartitions}
\alias{inferredpartitions}

\title{Identify inferred partitions in mapping QTL to a phylogenetic tree}

\description{
  Identify the inferred partitions for a chromosome from the results of scanPhyloQTL.
}

\usage{
inferredpartitions(output, chr, lodthreshold, probthreshold=0.9)
}
\arguments{
 \item{output}{An object output by the function
   \code{\link{scanPhyloQTL}}.}
 \item{chr}{A character string indicating the chromosome to
   consider. (It can also be a number, but it's then converted to a
   character string.)}
 \item{lodthreshold}{LOD threshold; if maximum LOD score is less than
   this, the null model is considered.}
 \item{probthreshold}{Threshold on posterior probabilities.  See Details
   below.} 
}


\details{
  We consider a single chromosome, and take the maximum LOD score for
  each partition on that chromosome.  The presence of a QTL is inferred
  if at least one partition has LOD score greater than
  \code{lodthreshold}.  In this case, we then convert the LOD scores for
  the partitions to approximate posterior probabilities by taking
  \eqn{10^{LOD}}{10^LOD} and then rescaling them to sum to 1.
  These are sorted from largest to smallest, and we 
  then take as the inferred partitions the smallest set whose posterior
  probabilities cumulatively add up to at least \code{probthreshold}.
}

\value{
  A vector of character strings.  If the null model (no QTL) is
  inferred, the output is \code{"null"}.  Otherwise, it is the set of
  inferred partitions.
}

\references{
  Broman, K. W., Kim, S., An\'e, C. and Payseur, B. A.  Mapping
  quantitative trait loci to a phylogenetic tree.  In preparation.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
\dontshow{set.seed(66787403)}
# example map; drop X chromosome
data(map10)           
map10 <- map10[1:19]
\dontshow{map10 <- map10[1:4]}

# simulate data
x <- simPhyloQTL(4, partition="AB|CD", crosses=c("AB", "AC", "AD"),
                 map=map10, n.ind=150,
                 model=c(1, 50, 0.5, 0))

# run calc.genoprob on each cross
\dontrun{x <- lapply(x, calc.genoprob, step=2)
}\dontshow{x <- lapply(x, calc.genoprob, step=0)}

# scan genome, at each position trying all possible partitions
out <- scanPhyloQTL(x, method="hk")

# inferred partitions
inferredpartitions(out, chr=3, lodthreshold=3)

# inferred partitions with prob'y threshold = 0.95
inferredpartitions(out, chr=3, lodthreshold=3, probthreshold=0.95)
}

\seealso{ \code{\link{scanPhyloQTL}}, \code{\link{plot.scanPhyloQTL}}, 
  \code{\link{summary.scanPhyloQTL}}, \code{\link{max.scanPhyloQTL}},
  \code{\link{simPhyloQTL}} }

\keyword{print}
