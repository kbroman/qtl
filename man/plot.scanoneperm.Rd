\name{plot.scanoneperm}
\alias{plot.scanoneperm}

\title{Plot permutation results for a single-QTL genome scan}

\description{
  Plot a histogram of the permutation results from a single-QTL genome
  scan.
}

\usage{
\method{plot}{scanoneperm}(x, lodcolumn=1, \dots)
}

\arguments{
  \item{x}{An object of class \code{"scanoneperm"}, as output by
    \code{\link{scanone}} when \code{n.perm} is specified.}
  \item{lodcolumn}{This indicates the LOD score column to plot.  This
    should be a single number between 1 and the number of LOD columns in
    the object input.}
  \item{\dots}{Passed to the function \code{\link[graphics]{hist}} when it
    is called.}
}
\value{None.}

\details{
  The function plots a histogram of the permutation results obtained by
  \code{\link{scanone}} when \code{n.perm} is specified. If
  separate permutations were performed for the autosomes and the X
  chromosome (using \code{perm.Xsp=TRUE}), separate histograms are given.
}

\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc)

operm <- scanone(fake.bc, method="hk", n.perm=100)
plot(operm)
}

\seealso{ \code{\link{scanone}},
  \code{\link{summary.scanoneperm}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
