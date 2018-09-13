\name{calc.penalties}
\alias{calc.penalties}

\title{Calculate LOD penalties}

\description{
  Derive penalties for the penalized LOD scores (used by
  \code{\link{stepwiseqtl}}) on the basis of permutation results
  from a two-dimensional, two-QTL scan (obtained by \code{\link{scantwo}}).
}

\usage{calc.penalties(perms, alpha=0.05, lodcolumn)}

\arguments{
 \item{perms}{Permutation results from \code{\link{scantwo}}.}
 \item{alpha}{Significance level.}
 \item{lodcolumn}{If the scantwo permutation results contain LOD scores for
    multiple phenotypes, this argument indicates which to use in the
    summary.  This may be a vector.  If missing, penalties for all
 phenotypes are calculated.}
}

\details{
  Thresholds derived from \code{\link{scantwo}} permutations (that
  is, for a two-dimensional, two-QTL genome scan) are used to calculate
  penalties on main effects and interactions.

  The main effect penalty is the 1-\code{alpha} quantile of the null
  distribution of the genome-wide maximum LOD score from a single-QTL
  genome scan (as with \code{\link{scanone}}).

  The "heavy" interaction penalty is the 1-\code{alpha} quantile of
  the null distribution of the maximum interaction LOD score (that is,
  the \eqn{\log_{10}}{log10} likelihood ratio comparing the best model
  with two interacting QTL to the best model with two additive QTL) from
  a two-dimensional, two-QTL genome scan (as with
  \code{\link{scantwo}}).

  The "light" interaction penality is the difference between the
  \code{"fv1"} threshold from the \code{\link{scantwo}}
  permutations (that is, the 1-\code{alpha} quantile of the LOD score
  comparing the best model with two interacting QTL to the best
  single-QTL model) and the main effect penalty.

  If the permutations results were obtained with \code{perm.Xsp=TRUE},
  to give X-chr-specific results, six penalties are calculated: main
  effect for autosomes, main effect for X chr, heavy penalty on A:A
  interactions, light penalty on A:A interactions, penalty on A:X
  interactions, and penalty on X:X interactions.
}

\value{
  Vector of three values indicating the penalty on main
  effects and heavy and light penalties on interactions,
  or a matrix of such results, with each row corresponding to a
  different phenotype.

  If the input permutations are X-chromosome-specific, the result has
  six values: main effect for autosomes, main effect for X chr, heavy
  penalty on A:A interactions, light penalty on A:A interactions,
  penalty on A:X interactions, and penalty on X:X interactions.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
fake.f2 <- calc.genoprob(fake.f2, step=5)
out.2dim <- scantwo(fake.f2, method="hk")

# permutations
\dontshow{permo.2dim <- scantwo(fake.f2, method="hk", n.perm=2)}
\dontrun{permo.2dim <- scantwo(fake.f2, method="hk", n.perm=1000)}
summary(permo.2dim, alpha=0.05)

# penalties
calc.penalties(permo.2dim)
}

\references{
  Manichaikul, A., Moon, J. Y., Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}, Yandell, B. S. and Broman,
  K. W. (2009) A model selection approach for the identification of
  quantitative trait loci in experimental crosses, allowing
  epistasis. \emph{Genetics}, \bold{181}, 1077--1086.
}

\seealso{ \code{\link{scantwo}}, \code{\link{stepwiseqtl}} }

\keyword{utilities}
