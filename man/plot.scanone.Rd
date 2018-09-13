\name{plot.scanone}
\alias{plot.scanone}

\title{Plot LOD curves}

\description{
  Plot the LOD curve for a genome scan with a single-QTL model (the
  output of \code{\link{scanone}}).
}

\usage{
\method{plot}{scanone}(x, x2, x3, chr, lodcolumn=1, incl.markers=TRUE,
     xlim, ylim, lty=1, col=c("black","blue","red"), lwd=2,
     add=FALSE, gap=25, mtick = c("line", "triangle"),
     show.marker.names=FALSE, alternate.chrid=FALSE,
     bandcol=NULL, type="l", cex=1, pch=1, bg="transparent",
     bgrect=NULL, \dots)
}

\arguments{
  \item{x}{An object of class \code{"scanone"}, as output by
    \code{\link{scanone}}.}
  \item{x2}{Optional second \code{scanone} object.}
  \item{x3}{Optional third \code{scanone} object.}
  \item{chr}{Optional vector indicating the chromosomes to plot.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
  \item{lodcolumn}{An integer, or vector of 3 integers, indicating which
    of the LOD score columns should be plotted (generally this is 1).}
  \item{incl.markers}{Indicate whether to plot line segments at the
  marker locations.}
  \item{xlim}{Limits for x-axis (optional).}
  \item{ylim}{Limits for y-axis (optional).}
  \item{lty}{Line types; a vector of length 1 or 3.}
  \item{col}{Line colors; a vector of length 1 or 3.}
  \item{lwd}{Line widths; a vector of length 1 or 3.}
  \item{add}{If TRUE, add to a current plot.}
  \item{gap}{Gap separating chromosomes (in cM).}
  \item{mtick}{Tick mark type for markers (line segments or
    upward-pointing triangels).}
  \item{show.marker.names}{If TRUE, show the marker names along the x axis.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
     axis labels, so that they may be more easily distinguished.}
  \item{bandcol}{Optional color for alternating bands to indicate
  chromosomes.  If NULL (the default), no bands are plotted.  A good
  choice might be \code{bandcol="gray70"}.}
  \item{type}{Type of plot (see \code{\link[graphics]{plot}}): for
  example, \code{type="l"} for lines or \code{type="p"} for points only,
  may be of length 1 or 3.}
  \item{cex}{Point size expansion, for example if \code{type="p"} is used.
  May be of length 1 or 3.}
  \item{pch}{Point type, for example if \code{type="p"} is used. See \code{\link[graphics]{points}}.
  May be of length 1 or 3.}
  \item{bg}{Background color for points, for example if \code{type="p"}
  and \code{pch=21} are used. See \code{\link[graphics]{points}}.
  May be of length 1 or 3.}
  \item{bgrect}{Optional background color for the rectangular plotting region.}
  \item{\dots}{Passed to the function \code{\link[graphics]{plot}} when it
    is called.}
}
\value{None.}

\details{
  This function allows you to plot the results of up to three genome
  scans against one another.  Such objects must conform with each
  other.

  One may alternatively use the argument \code{add} to add the plot of
  an additional genome scan to the current figure, but some care is
  required: the same chromosomes should be selected, and the results
  must concern crosses with the same genetic maps.

  If a single \code{scanone} object containing multiple LOD score
  columns (for example, from different phenotypes) is input, up to three
  LOD curves may be plotted, by providing a vector in the argument
  \code{lodcolumn}.  If multiple \code{scanone} objects are input (via
  \code{x}, \code{x2} and \code{x3}), the LOD score columns to be
  plotted are chosen from the corresponding element of the
  \code{lodcolumn} argument.
}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(1,13))}
fake.f2 <- calc.genoprob(fake.f2,step=2.5)
out.mr <- scanone(fake.f2, method="mr")
out.em <- scanone(fake.f2, method="em")
plot(out.mr)
plot(out.mr, out.em, chr=c(1,13), lty=1, col=c("violetred","black"))
out.hk <- scanone(fake.f2, method="hk")
plot(out.hk, chr=c(1,13), add=TRUE, col="slateblue")

plot(out.hk, chr=13, show.marker.names=TRUE)

plot(out.hk, bandcol="gray70")

# plot points rather than lines
plot(out.hk, bandcol="gray70", type="p", cex=0.3, pch=21, bg="slateblue")
}

\seealso{ \code{\link{scanone}},
  \code{\link{summary.scanone}}, \code{\link[graphics]{par}},
  \code{\link[grDevices]{colors}}, \code{\link{add.threshold}}, \code{\link{xaxisloc.scanone}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
