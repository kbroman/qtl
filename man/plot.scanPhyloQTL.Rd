\name{plot.scanPhyloQTL}
\alias{plot.scanPhyloQTL}

\title{Plot LOD curves from single-QTL scan to map QTL to a phylogenetic tree}

\description{
  Plot the LOD curves for each partition for a genome scan with a single
  diallelic QTL (the
  output of \code{\link{scanPhyloQTL}}).
}

\usage{
\method{plot}{scanPhyloQTL}(x, chr, incl.markers=TRUE,
     col, xlim, ylim, lwd=2, gap=25, mtick=c("line", "triangle"),
     show.marker.names=FALSE, alternate.chrid=FALSE, legend=TRUE, \dots)
}

\arguments{
  \item{x}{An object of class \code{"scanPhyloQTL"}, as output by
    \code{\link{scanPhyloQTL}}.}
  \item{chr}{Optional vector indicating the chromosomes to plot.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
  \item{incl.markers}{Indicate whether to plot line segments at the
  marker locations.}
  \item{col}{Optional vector of colors to use for each partition.}
  \item{xlim}{Limits for x-axis (optional).}
  \item{ylim}{Limits for y-axis (optional).}
  \item{lwd}{Line width.}
  \item{gap}{Gap separating chromosomes (in cM).}
  \item{mtick}{Tick mark type for markers (line segments or
    upward-pointing triangels).}
  \item{show.marker.names}{If TRUE, show the marker names along the x axis.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
     axis labels, so that they may be more easily distinguished.}
  \item{legend}{Indicates whether to include a legend in the plot.}
  \item{\dots}{Passed to the function \code{\link{plot.scanone}} when it
    is called.}
}
\value{None.}

%\details{
%}

\examples{
\dontrun{
# example map; drop X chromosome
data(map10)
map10 <- map10[1:19]

# simulate data
x <- simPhyloQTL(4, partition="AB|CD", crosses=c("AB", "AC", "AD"),
                 map=map10, n.ind=150,
                 model=c(1, 50, 0.5, 0))

# run calc.genoprob on each cross
x <- lapply(x, calc.genoprob, step=2)

# scan genome, at each position trying all possible partitions
out <- scanPhyloQTL(x, method="hk")

# maximum peak
max(out, format="lod")

# approximate posterior probabilities at peak
max(out, format="postprob")

# all peaks above a threshold for LOD(best) - LOD(2nd best)
summary(out, threshold=1, format="lod")

# all peaks above a threshold for LOD(best), showing approx post'r prob
summary(out, format="postprob", threshold=3)

# plot of results
plot(out)
}
}

\seealso{ \code{\link{scanPhyloQTL}}, \code{\link{max.scanPhyloQTL}},
  \code{\link{summary.scanPhyloQTL}}, \code{\link{plot.scanone}},
  \code{\link{inferredpartitions}},
  \code{\link{simPhyloQTL}},
  \code{\link[graphics]{par}}, \code{\link[grDevices]{colors}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\references{
  Broman, K. W., Kim, S., An\'e, C. and Payseur, B. A.  Mapping
  quantitative trait loci to a phylogenetic tree.  In preparation.
}

\keyword{hplot}
