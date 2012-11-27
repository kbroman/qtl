\name{xaxisloc.scanone}
\alias{xaxisloc.scanone}

\title{Get x-axis locations in scanone plot}

\description{
  Get x-axis locations for given cM positions on given chromosomes in a
  plot from \code{\link{plot.scanone}})
}

\usage{
xaxisloc.scanone(out, thechr, thepos, chr, gap=25)
}			    

\arguments{
  \item{out}{An object of class \code{"scanone"}, as output by
    \code{\link{scanone}}.  This must be identical to what was used
  in the call to \code{\link{plot.scanone}}.}
  \item{thechr}{Chromosome IDs at which x-axis locations are to be determined.}
  \item{thepos}{Chromosome positions at which x-axis locations are to be
    determined.}
  \item{chr}{Optional vector specifying which chromosomes were plotted.
  This must be identical to what was used in the call to
  \code{\link{plot.scanone}}.} 
  \item{gap}{Gap separating chromosomes (in cM).  This must be identical
  to what was used in the call to \code{\link{plot.scanone}}.}
}

\value{A numeric vector of x-axis locations.}

\details{
  This function allows you to identify the x-axis locations in a plot of 
  genome scan results, produced by
  \code{\link{plot.scanone}}.  This is useful for adding
  annotations, such as text or arrows.

  The arguments \code{out}, \code{chr}, and \code{gap} must match what
  was used in the call to \code{\link{plot.scanone}}.

  The arguments \code{thechr} and \code{thepos} indicate the genomic
  positions for which x-axis locations are desired.  If they both have
  length > 1, they must have the same length.  If one has length > 1 and
  one has length 1, the one with length 1 is expanded to match.
}


\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, c(1,4,6,15))}
hyper <- calc.genoprob(hyper)
out <- scanone(hyper, method="hk")
plot(out, chr=c(1, 4, 6, 15))

# add arrow and text to indicate peak LOD score
mxout <- max(out)
x <- xaxisloc.scanone(out, mxout$chr, mxout$pos, chr=c(1,4,6,15))
arrows(x+30, mxout$lod, x+5, mxout$lod, len=0.1, col="blue")
text(x+35, mxout$lod, "the peak", col="blue", adj=c(0, 0.5))
}

\seealso{ 
  \code{\link{plot.scanone}}, 
  \code{\link{add.threshold}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
