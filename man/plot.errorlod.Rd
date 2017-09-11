\name{plotErrorlod}
\alias{plotErrorlod}

\title{Plot grid of error LOD values}

\description{
  Plot a grid of the LOD scores indicating which genotypes are
  likely to be in error.
}

\usage{
plotErrorlod(x, chr, ind, breaks=c(-Inf,2,3,4.5,Inf),
              col=c("white","gray85","hotpink","purple3"),
              alternate.chrid=FALSE, \dots)
}

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link{read.cross}} for details.}
  \item{chr}{Optional vector indicating the chromosomes to be drawn in
  the plot. This should be a vector of character strings referring to
  chromosomes by name; numeric values are converted to strings.  Refer
  to chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
  \item{ind}{Indicates the individuals for which the error LOD scores
    should be plotted (passed to \code{\link{subset.cross}}).}
  \item{breaks}{A set of breakpoints for the colors; must give one more
    breakpoint than color.  Intervals are open on the left and closed on
    the right, except for the lowest interval.}
  \item{col}{A vector of colors to appear in the image.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
     axis labels, so that they may be more easily distinguished.}
  \item{\dots}{Ignored at this point.}
}

\details{
  Uses \code{\link[graphics]{image}} to plot a grid with different shades
  of pixels to indicate which genotypes are likely to be in error.

  Darker pixels have higher error LOD scores:
  \eqn{LOD \le 2} in white;
  \eqn{2 < LOD \le 3}{2 < LOD <= 3} in gray;
  \eqn{3 < LOD \le 4.5}{3 < LOD <= 4.5} in pink;
  \eqn{LOD > 4.5} in purple.
}

\value{None.}

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper,chr=c(1,19))}
# Calculate error LOD scores
hyper <- calc.errorlod(hyper,error.prob=0.01)

# plot the error LOD scores; print those above a specified cutoff
plotErrorlod(hyper)
plotErrorlod(hyper,chr=1)
}

\seealso{ \code{\link{calc.errorlod}},
  \code{\link{top.errorlod}}, \code{\link[graphics]{image}},
  \code{\link{subset.cross}}, \code{\link{plotGeno}} }

\references{
  Lincoln, S. E. and Lander, E. S. (1992) Systematic detection of
  errors in genetic linkage data.  \emph{Genomics} \bold{14}, 604--610.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
