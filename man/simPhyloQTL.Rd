\name{simPhyloQTL}
\alias{simPhyloQTL}

\title{Simulate a set of intercrosses for a single diallelic QTL}

\description{
  Simulate a set of intercrosses with a single diallelic QTL.
}

\usage{
simPhyloQTL(n.taxa=3, partition, crosses, map, n.ind=100, model,
            error.prob=0, missing.prob=0, partial.missing.prob=0,
            keep.qtlgeno=FALSE, keep.errorind=TRUE, m=0, p=0,
  	    map.function=c("haldane","kosambi","c-f","morgan"))
}

\arguments{
 \item{n.taxa}{Number of taxa (i.e., strains).}
 \item{partition}{A vector of character strings of the form "AB|CD" or "A|BCD"
 indicating, for each QTL, which taxa have which allele.  If missing,
 simulate under the null hypothesis of no QTL.}
 \item{crosses}{A vector of character strings indicating the crosses to
 do (for the form "AB", "AC", etc.).  These will be sorted and then only
 unique ones used.  If missing, all crosses will be simulated.}
 \item{map}{A list whose components are vectors containing the marker
   locations on each of the chromosomes.}
 \item{n.ind}{The number of individuals in each cross.  If length 1, all
 crosses will have the same number of individuals; otherwise the length
 should be the same as \code{crosses}.}
 \item{model}{A matrix where each row corresponds to a
 different QTL, and gives the chromosome number, cM position and
 effects of the QTL (assumed to be the same in each cross in which the
 QTL is segregating).}
 \item{error.prob}{The genotyping error rate.}
 \item{missing.prob}{The rate of missing genotypes.}
 \item{partial.missing.prob}{When simulating an intercross or 4-way
   cross, this gives the rate at which markers will be incompletely
   informative (i.e., dominant or recessive).}
 \item{keep.qtlgeno}{If TRUE, genotypes for the simulated QTLs will be
   included in the output.}
 \item{keep.errorind}{If TRUE, and if \code{error.prob > 0}, the
   identity of genotyping errors will be included in the output.}
 \item{m}{Interference parameter; a non-negative integer.  0 corresponds
   to no interference.}
 \item{p}{Probability that a chiasma comes from the no-interference
   mechanism} 
 \item{map.function}{Indicates whether to use the Haldane, Kosambi, 
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
}

\details{
  Meiosis is assumed to follow the Stahl model for crossover
  interference (see the references, below), of which the no interference
  model and the chi-square model are special cases.  Chiasmata on the
  four-strand bundle are a superposition of chiasmata from two different
  mechanisms.  With probability \code{p}, they arise by a mechanism
  exhibiting no interference; the remainder come from a chi-square model
  with inteference parameter \code{m}.  Note that \code{m=0} corresponds
  to no interference, and with \code{p=0}, one gets a pure chi-square
  model.
  
  QTLs are assumed to act additively, and the residual phenotypic
  variation is assumed to be normally distributed with variance 1.

  The effect of a QTL is a pair of numbers,
  (\eqn{a,d}), where \eqn{a} is the additive effect (half the difference
  between the homozygotes) and \eqn{d} is the dominance deviation (the
  difference between the heterozygote and the midpoint between the
  homozygotes).

}  

\value{
  A list with each component being an object of class \code{cross}.  See \code{\link{read.cross}} for
  details.  The names (e.g. "AB", "AC", "BC") indicate the crosses.

  If \code{keep.qtlgeno} is TRUE, each cross object will contain a
  component \code{qtlgeno} which is a matrix containing the QTL
  genotypes (with complete data and no errors), coded as in the genotype
  data.

  If \code{keep.errorind} is TRUE and errors were simulated, each
  component of \code{geno} in each cross will each contain a matrix \code{errors},
  with 1's indicating simulated genotyping errors.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\seealso{ \code{\link{scanPhyloQTL}},   \code{\link{inferredpartitions}}, 
  \code{\link{summary.scanPhyloQTL}}, \code{\link{max.scanPhyloQTL}},
  \code{\link{plot.scanPhyloQTL}}, 
  \code{\link{sim.cross}}, \code{\link{read.cross}} }

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
                    
\references{
  Broman, K. W., Kim, S., An\'e, C. and Payseur, B. A.  Mapping
  quantitative trait loci to a phylogenetic tree.  In preparation.
}
  
\keyword{datagen}
