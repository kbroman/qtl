\name{scantwopermhk}
\alias{scantwopermhk}

\title{Permutation test for 2d genome scan by Haley-Knott regression}

\description{
  Perform a permutation test with a two-dimensional genome scan with a
  two-QTL model, with possible allowance for additive covariates, by
  Haley-Knott regression.
}

\usage{
scantwopermhk(cross, chr, pheno.col=1,
              addcovar=NULL, weights=NULL, n.perm=1,
              batchsize=1000,
              perm.strata=NULL, perm.Xsp=NULL,
              verbose=FALSE, assumeCondIndep=FALSE)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes for which LOD
   scores should be calculated.  This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This should be a single value (numeric index
   or character string for a phenotype name), but it may also be a
   vector of numeric values with length equal to the number
   of individuals in the cross, in which case it is taken to be a vector
   of individuals' phenotypes.}
 \item{addcovar}{Additive covariates.}
 \item{weights}{Optional weights of individuals.  Should be either NULL
   or a vector of length n.ind containing positive weights. Used only
   in the case \code{model="normal"}.}
 \item{n.perm}{Number of permutation replicates.}
 \item{batchsize}{If \code{n.perm > batchsize}, permutations will be
   run in batches of no more than \code{batchsize} permutations.}
 \item{perm.strata}{Used to perform a
   stratified permutation test.  This should be a vector with the same
   number of individuals as in the cross data.  Unique values indicate
   the individual strata, and permutations will be performed within the
   strata.}
 \item{perm.Xsp}{If TRUE, run separate permutations for A:A, A:X, and
 X:X. In this case, \code{n.perm} refers to the number of permutations
 for the A:A part; more permutations are used for the A:X and X:X parts,
 as estimates of quantiles farther out into the tails are needed.}
 \item{verbose}{If TRUE, display information about the progress of
   calculations.}
 \item{assumeCondIndep}{If TRUE, assume conditional independence of QTL
   genotypes given marker genotypes.  This is an approximation, but it
   may speed things up.}
}

\details{
  This is a scaled-back version of the permutation test provided by
  \code{\link{scantwo}}: only for a normal model with Haley-Knott
  regression, and not allowing interactive covariates.

  This is an attempt to speed things up and attentuate the memory usage
  problems in \code{\link{scantwo}}.

  In the case of \code{perm.Xsp=TRUE} (X-chr-specific thresholds), we
  use a stratified permutation test, stratified by sex and
  cross-direction.
}

\value{
  A list with six
  different LOD scores from each of the permutation replicates.
  First, the maximum LOD score for the full model (two QTLs plus an
  interaction).  Second, for each pair of
  chromosomes, we take the difference between the full LOD and the
  maximum single-QTL LOD for those two chromosomes, and then maximize
  this across chromosome pairs.  Third, for each pair of chromosomes we
  take the difference between the maximum full LOD and the maximum
  additive LOD, and then maximize this across chromosome pairs.  Fourth,
  the maximum LOD score for the additive QTL model.  Fifth, for each
  pair of chromosomes, we take the difference between the additive LOD
  and the maximum single-QTL LOD for those two chromosomes, and then
  maximize this across chromosome pairs.  Finally, the maximum
  single-QTL LOD score (that is, from a single-QTL scan).  The latter is
  not used in \code{\link{summary.scantwoperm}}, but does get
  calculated at each permutation, so we include it for the sake of
  completeness.

  If \code{perm.Xsp=TRUE}, this is a list of lists, for the A:A, A:X,
  and X:X sections, each being a list as described above.
}

\references{
  Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for
  quantitative trait mapping.  \emph{Genetics} \bold{138}, 963--971.

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.
}

\author{Karl W Broman, \email{broman@wisc.edu}; Hao Wu}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
fake.f2 <- calc.genoprob(fake.f2, step=5)
operm <- scantwopermhk(fake.f2, n.perm=2)
summary(operm, alpha=0.05)
}

\seealso{ \code{\link{scantwo}}, \code{\link{plot.scantwoperm}},
  \code{\link{summary.scantwoperm}},
  \code{\link{c.scantwoperm}} }

\keyword{models}
