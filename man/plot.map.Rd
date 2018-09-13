\name{plotMap}
\alias{plotMap}
\alias{plot.map}

\title{Plot genetic map}

\description{
  Plot genetic map of marker locations for all chromosomes.
}

\usage{
\method{plot}{map}(x, map2, chr, horizontal=FALSE, shift=TRUE,
     show.marker.names=FALSE, alternate.chrid=FALSE, \dots)
plotMap(x, map2, chr, horizontal=FALSE, shift=TRUE,
     show.marker.names=FALSE, alternate.chrid=FALSE, \dots)
}

\arguments{
  \item{x}{A list whose components are vectors of marker locations.  A
    \code{cross} object may be given instead, in which case the genetic
    map it contains is used.}
  \item{map2}{An optional second genetic map with the same number (and
    names) of chromosomes.  As with the first argument, a
    \code{cross} object may be given instead.  If this argument is
    given, a comparison of the two genetic maps is plotted.}
  \item{chr}{Optional vector indicating the chromosomes to plot.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
  \item{horizontal}{Specifies whether the chromosomes should be
    plotted horizontally.}
  \item{shift}{If TRUE, shift the first marker on each chromosome to be
    at 0 cM.}
  \item{show.marker.names}{If TRUE, marker names are included.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
     axis labels, so that they may be more easily distinguished.}
  \item{\dots}{Passed to \code{\link[graphics]{plot}}.}
}

\value{None.}

\details{
  Plots the genetic map for each chromosome, or a comparison of the
  genetic maps if two maps are given.

  For a comparison of two maps, the first
  map is on the left (or, if \code{horizontal=TRUE}, on the top).  Lines
  are drawn to connect markers.  Markers that exist in just one map and
  not the other are indicated by short line segments, on one side or the
  other, that are not connected across.

  For a sex-specific map, female and male maps are plotted against one
  another.  For two sex-specific maps, the two female maps are plotted
  against one another and the two male maps are plotted against one
  another.
}

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=18:19)}
plotMap(fake.bc)
plotMap(fake.bc,horizontal=TRUE)

newmap <- est.map(fake.bc)
plot(newmap)
plotMap(fake.bc, newmap)

plotMap(fake.bc, show.marker.names=TRUE)

}

\seealso{ \code{\link{est.map}}, \code{\link{plot.cross}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
