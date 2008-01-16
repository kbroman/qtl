\name{plot.scanoneperm}
\alias{plot.scanoneperm}

\title{Plot permutation results for a single-QTL genome scan}

\description{
  Plot a histogram of the permutation results from a single-QTL genome
  scan. 
}

\usage{
\method{plot}{scanoneperm}(x, \dots)
}			    

\arguments{
  \item{x}{An object of class \code{"scanoneperm"}, as output by
    \code{\link[qtl]{scanone}} when \code{n.perm} is specified.}
  \item{\dots}{Passed to the function \code{\link[graphics]{hist}} when it
    is called.}
}
\value{None.}

\details{
  The function plots a histogram of the permutation results obtained by
  \code{\link[qtl]{scanone}} when \code{n.perm} is specified. If
  separate permutations were performed for the autosomes and the X
  chromosome (using \code{perm.Xsp=TRUE}), separate histograms are given.
}

\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc)

operm <- scanone(fake.bc, method="hk", n.perm=100)
plot(operm)
}

\seealso{ \code{\link[qtl]{scanone}},
  \code{\link[qtl]{summary.scanoneperm}} } 

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
