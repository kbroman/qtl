\name{summary.scanoneperm}
\alias{summary.scanoneperm}

\title{LOD thresholds from scanone permutation results}

\description{
  Print the estimated genome-wide LOD thresholds on the basis of
  permutation results from \code{\link{scanone}} (with
  \code{n.perm} > 0).
}

\usage{
\method{summary}{scanoneperm}(object, alpha=c(0.05, 0.10),
                    controlAcrossCol=FALSE, \dots)
}

\arguments{
 \item{object}{Output from the function \code{\link{scanone}}
   with \code{n.perm} > 0.}
 \item{alpha}{Genome-wide significance levels.}
 \item{controlAcrossCol}{If TRUE, control error rate not just across
   the genome but also across the columns of LOD scores.}
 \item{\dots}{Ignored at this point.}
}

\details{
  If there were autosomal data only or \code{\link{scanone}} was
  run with \code{perm.Xsp=FALSE}, genome-wide LOD thresholds are given;
  these are the 1-\eqn{\alpha}{alpha} quantiles of the genome-wide maximum LOD
  scores from the permutations.

  If there were autosomal and X chromosome data and
  \code{\link{scanone}} was run with \code{perm.Xsp=TRUE},
  autosome- and X-chromsome-specific LOD thresholds are given, by the
  method described in Broman et al. (2006).  Let \eqn{L_A} and
  \eqn{L_X} be total the genetic lengths of the autosomes and X
  chromosome, respectively, and let \eqn{L_T = L_A + L_X}
  Then in place of \eqn{\alpha}{alpha}, we use
  \deqn{\alpha_A = 1 - (1-\alpha)^{L_A/L_T}}{alpha_A = 1 - (1 -
    alpha)^(L_A/L_T)}
  as the significance level for the autosomes and
  \deqn{\alpha_X = 1 - (1-\alpha)^{L_X/L_T}}{alpha_x = 1 - (1 -
    alpha)^(LX/LT)}
  as the significance level for the X chromosome.  The result is a list
  with two matrices, one for the autosomes and one for the X chromosome.

  If \code{controlAcrossCol=TRUE}, we use a trick to control the error
  rate not just across the genome but also across the LOD score
  columns.  Namely, we convert each column of permutation results to
  ranks, and then for each permutation replicate we find the maximum
  rank across the columns.  We then find the appropriate quantile of the
  maximized ranks, and then backtrack to the corresponding LOD score
  within each of the columns.
}

\value{
  An object of class \code{summary.scanoneperm}, to be printed by
  \code{print.summary.scanoneperm}.  If there were X chromosome data and
  \code{scanone} was run with \code{perm.Xsp=TRUE}, there are two
  matrices in the results, for the autosome and X-chromosome LOD
  thresholds.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"))}
fake.f2 <- calc.genoprob(fake.f2, step=2.5)

operm1 <- scanone(fake.f2, n.perm=100, method="hk")
summary(operm1)

operm2 <- scanone(fake.f2, n.perm=100, method="hk", perm.Xsp=TRUE)
summary(operm2)

# Add noise column
fake.f2$pheno$noise <- rnorm(nind(fake.f2))
operm3 <- scanone(fake.f2, pheno.col=c("phenotype", "noise"), n.perm=10, method="hk")
summary(operm3)
summary(operm3, controlAcrossCol=TRUE, alpha=c(0.05, 0.36))
}

\seealso{ \code{\link{scanone}},
  \code{\link{summary.scanone}},
  \code{\link{plot.scanoneperm}} }

\references{
  Broman, K. W., Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}, Owens, S. E., Manichaikul, A.,
  Southard-Smith, E. M. and Churchill G. A. (2006)  The X chromosome in
  quantitative trait locus mapping.  \emph{Genetics}, \bold{174}, 2151--2158.

  Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for
  quantitative trait mapping.  \emph{Genetics} \bold{138}, 963--971.
}


\keyword{print}
