\name{plot.scantwoperm}
\alias{plot.scantwoperm}

\title{Plot permutation results for a 2d, 2-QTL genome scan}

\description{
  Plot a histogram of the permutation results from a two-dimensional,
  two-QTL genome scan.
}

\usage{
\method{plot}{scantwoperm}(x, lodcolumn=1, include_rug=TRUE, \dots)
}

\arguments{
  \item{x}{An object of class \code{"scantwoperm"}, as output by
    \code{\link{scantwo}} when \code{n.perm} is specified.}
  \item{lodcolumn}{This indicates the LOD score column to plot.  This
    should be a single number between 1 and the number of LOD columns in
    the object input.}
  \item{include_rug}{If TRUE, include a call to \code{\link[graphics]{rug}}.}
  \item{\dots}{Passed to the function \code{\link[graphics]{hist}} when it
    is called.}
}
\value{None.}

\details{
  The function plots a histogram of the permutation results obtained by
  \code{\link{scantwo}} when \code{n.perm} is specified. Separate
  histograms are provided for the five LOD scores, \code{full},
  \code{fv1}, \code{int}, \code{add}, and \code{av1}.
}

\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc)

operm2 <- scantwo(fake.bc, method="hk", n.perm=10)
plot(operm2)
}

\seealso{ \code{\link{scantwo}},
  \code{\link{summary.scantwoperm}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
