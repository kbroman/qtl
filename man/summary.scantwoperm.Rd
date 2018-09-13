\name{summary.scantwoperm}
\alias{summary.scantwoperm}

\title{LOD thresholds from scantwo permutation results}

\description{
  Print the estimated genome-wide LOD thresholds on the basis of
  permutation results from \code{\link{scantwo}} (with
  \code{n.perm} > 0).
}

\usage{
\method{summary}{scantwoperm}(object, alpha=c(0.05, 0.10), \dots)
}

\arguments{
 \item{object}{Output from the function \code{\link{scantwo}}
   with \code{n.perm} > 0.}
 \item{alpha}{Genome-wide significance levels.}
 \item{\dots}{Ignored at this point.}
}

\details{
  We take the \eqn{1-\alpha}{1-alpha} quantiles of the individual LOD
  scores.

  In the case of X-chr-specific permutations, we use the combined length
  of the autosomes, \eqn{L_A}, and the length of the X chromosome,
  \eqn{L_X}, and calculate the area of the A:A, A:X, and X:X regions as
  \eqn{L_A^2/2}, \eqn{L_A L_X}, and \eqn{L_X^2/2}, and then use the
  nominal significance levels of \eqn{1 - (1-\alpha)^p}{1 - (1-alpha)^p},
  where \eqn{p} is the proportional area for that region.
}

\value{
  An object of class \code{summary.scantwoperm}, to be printed by
  \code{print.summary.scantwoperm}.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"), ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=0)

\dontrun{operm <- scantwo(fake.f2, n.perm=100, method="hk")}\dontshow{
operm <- scantwo(fake.f2, n.perm=3, method="hk")}
summary(operm)
}

\seealso{ \code{\link{scantwo}},
  \code{\link{summary.scantwo}},
  \code{\link{plot.scantwoperm}} }

\references{
  Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for
  quantitative trait mapping.  \emph{Genetics} \bold{138}, 963--971.
}


\keyword{print}
