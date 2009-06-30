\name{add.threshold}
\alias{add.threshold}

\title{Add significance threshold to plot}

\description{
  Add a significance threshold to a plot created by 
  \code{\link[qtl]{plot.scanone}}), using the permutation results.
}

\usage{
add.threshold(out, chr, perms, alpha=0.05, lodcolumn=1, gap=25, \dots)
}			    

\arguments{
  \item{out}{An object of class \code{"scanone"}, as output by
    \code{\link[qtl]{scanone}}.  This must be identical to what was used
  in the call to \code{\link[qtl]{plot.scanone}}.}
  \item{chr}{Optional vector specifying which chromosomes to plot. If a
    selected subset of chromosomes were plotted, they must be specified
    here.} 
  \item{perms}{Permutation results from \code{\link[qtl]{scanone}}, used
    to calculate the significance threshold.}
  \item{alpha}{Significance level of the threshold.}
  \item{lodcolumn}{An integer indicating which
    of column in the permutation results should be used.}
  \item{gap}{Gap separating chromosomes (in cM).  This must be identical
  to what was used in the call to \code{\link[qtl]{plot.scanone}}.}
  \item{\dots}{Passed to the function \code{\link[graphics]{abline}} when it
    is called.}
}
\value{None.}

\details{
  This function allows you to add a horizontal line at the significance
  threshold to genome scan results plotted by
  \code{\link[qtl]{plot.scanone}}.

  The arguments \code{out}, \code{chr}, and \code{gap} must match what
  was used in the call to \code{\link[qtl]{plot.scanone}}.

  The argument \code{perms} must be specified.  If X-chromosome-specific
  permutations were performed (via the argument \code{perm.Xsp} in the
  call to \code{\link[qtl]{scanone}}), separate thresholds will be
  plotted for the autosomes and the X chromosome.  These are calculated
  via the \code{\link[qtl]{summary.scanoneperm}} function.
}

\examples{
data(hyper)
hyper <- calc.genoprob(hyper)
out <- scanone(hyper, method="hk")
operm <- scanone(hyper, method="hk", n.perm=100, perm.Xsp=TRUE)

plot(out, chr=c(1,4,6,15,"X"))
add.threshold(out, chr=c(1,4,6,15,"X"), perms=operm, alpha=0.05)
add.threshold(out, chr=c(1,4,6,15,"X"), perms=operm, alpha=0.1,
              col="green", lty=2)
}

\seealso{ \code{\link[qtl]{scanone}},
  \code{\link[qtl]{plot.scanone}}, 
  \code{\link[qtl]{summary.scanoneperm}}, \code{\link[qtl]{xaxisloc.scanone}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
