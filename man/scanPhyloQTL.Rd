\name{scanPhyloQTL}
\alias{scanPhyloQTL}

\title{Single-QTL genome scan to map QTL to a phylogenetic tree}

\description{
  Jointly consider multiple intercrosses with a single diallelic QTL
  model, considering all possible partitions of the strains into the two
  QTL allele groups.
}

\usage{
scanPhyloQTL(crosses, partitions, chr, pheno.col=1,
             model=c("normal", "binary"), method=c("em", "imp", "hk"),
             addcovar, maxit=4000, tol=0.0001, useAllCrosses=TRUE,
             verbose=FALSE)
}
\arguments{
 \item{crosses}{A list with each component being an intercross, as an object of class
   \code{cross} (see \code{\link{read.cross}} for details).  The names
   (of the form "AB") indicate the strains in the cross.}
 \item{partitions}{A vector of character strings of the form "AB|CD" or
   "A|BCD" indicating the set of paritions of the strains into two
   allele groups.  If missing, all partitions should be considered.}
 \item{chr}{Optional vector indicating the chromosomes for which LOD
   scores should be calculated.  This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This can be a vector of integers; for methods
   \code{"hk"} and \code{"imp"} this can be considerably faster than doing
   them one at a time.  One may also give a character strings matching
   the phenotype names.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{model}{The phenotype model: the usual normal model or a model for
   binary traits}
 \item{method}{Indicates whether to use the EM algorithm, 
   imputation, or Haley-Knott regression.}
 \item{addcovar}{Optional set of additive covariates to include in the
   analysis, as a list with the same length as \code{crosses}.  They
   must be numeric vectors or matrices, as for \code{\link{scanone}}.}
 \item{maxit}{Maximum number of iterations for method \code{"em"}.}
 \item{tol}{Tolerance value for determining convergence for method
   \code{"em"}.}
 \item{useAllCrosses}{If TRUE, use all crosses in the analysis of all
   partitions, with crosses not segregating the QTL included in the
   estimation of the residual variance.}
 \item{verbose}{If TRUE, print information about progress.}
}


\details{
  The aim is to jointly consider multiple intercrosses to not just map
  QTL but to also, under the assumption of a single diallelic QTL,
  identify the set of strains with each QTL allele.

  For each partition (of the strains into two groups) that is under
  consideration, we pull out the set of crosses that are segregating the
  QTL, re-code the alleles, and combine the crosses into one large cross.
  Crosses not segregating the QTL are also used, though with no QTL
  effects.  

  Additive covariate indicators for the crosses are included in the
  analysis, to allow for the possibility that there are overall shifts
  in the phenotypes between crosses.
}

\value{
  A data frame, as for the output of \code{\link{scanone}}, though with
  LOD score columns for each partition that is considered.  The result
  is given class \code{"scanPhyloQTL"}.
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

# maximum peak
max(out, format="lod")

# approximate posterior probabilities at peak
max(out, format="postprob")

# all peaks above a threshold for LOD(best) - LOD(2nd best)
summary(out, threshold=1, format="lod")

# all peaks above a threshold for LOD(best), showing approx post'r prob
summary(out, format="postprob", threshold=3)

# plot results
plot(out)
}

\seealso{ \code{\link{plot.scanPhyloQTL}}, 
  \code{\link{summary.scanPhyloQTL}}, \code{\link{max.scanPhyloQTL}},
  \code{\link{inferredpartitions}}, 
  \code{\link{simPhyloQTL}} }

\keyword{models}
