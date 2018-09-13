\name{add.cim.covar}
\alias{add.cim.covar}

\title{Indicate marker covariates from composite interval mapping}

\description{
  Add dots at the locations of the selected marker covariates, for a
  plot of composite interval mapping results.
}

\usage{
add.cim.covar(cimresult, chr, gap=25, \dots)
}

\arguments{
 \item{cimresult}{Composite interval mapping results, as output from
   \code{\link{cim}}.} 
  \item{chr}{Optional vector specifying which chromosomes to plot. (The
    chromosomes must be specified by name.)  This should be identical to
  that used in the call to \code{\link{plot.scanone}}.}
  \item{gap}{Gap separating chromosomes (in cM). This should be identical to
    that used in the call to \code{\link{plot.scanone}}.}
  \item{\dots}{Additional plot arguments, passed to the function
    \code{\link[graphics]{points}}.} 
}


\details{
  One must first have used the function \code{\link{plot.scanone}}
  to plot the composite interval mapping results.

  The arguments \code{chr} and \code{gap} must be identical to the
  values used in the call to \code{\link{plot.scanone}}.

  Dots indicating the locations of the selected marker covariates are
  displayed on the x-axis.  (By default, solid red circles are plotted;
  this may be modified by specifying the graphics parameters \code{pch} and
  \code{col}.)
}

\value{
  A data frame indicating the marker covariates that were plotted.
}


\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
\dontrun{data(hyper)
hyper <- calc.genoprob(hyper, step=2.5)

out <- scanone(hyper)
out.cim <- cim(hyper, n.marcovar=3)
plot(out, out.cim, chr=c(1,4,6,15), col=c("blue", "red"))

add.cim.covar(out.cim, chr=c(1,4,6,15))}
}

\seealso{ \code{\link{cim}}, \code{\link{plot.scanone}} }

\keyword{hplot}
