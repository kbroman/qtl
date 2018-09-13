\name{plot.cross}
\alias{plot.cross}

\title{Plot various features of a cross object}

\description{
  Plots grid of the missing genotypes, genetic map, and histograms or
  barplots of phenotypes for the data from an experimental cross.
}

\usage{
\method{plot}{cross}(x, auto.layout=TRUE, pheno.col,
          alternate.chrid=TRUE, \dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{auto.layout}{If TRUE, \code{par(mfrow)} is set so that all plots
   fit within one figure.}
 \item{pheno.col}{Vector of numbers or character strings corresponding to
   phenotypes that should be plotted. If unspecified, all phenotypes are
   plotted.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
     axis labels, so that they may be more easily distinguished.}
 \item{\dots}{Ignored at this point.}
}

\value{None.}

\details{
  Calls \code{\link{plotMissing}}, \code{\link{plotMap}} and
  \code{\link{plotPheno}} to plot the missing genotypes, genetic
  map, and histograms or barplots of all phenotypes.

  If \code{auto.format=TRUE}, \code{par(mfrow)} is used with
  \code{ceiling(sqrt(n.phe+2))} rows and the minimum number of columns
  so that all plots fit on the plotting device.

  Numeric phenotypes are displayed as histograms or barplots by calling
  \code{\link{plotPheno}}.
}

\examples{
data(fake.bc)
plot(fake.bc)
}

\seealso{ \code{\link{plotMissing}}, \code{\link{plotMap}},
  \code{\link{plotPheno}} }


\author{Karl W Broman, \email{broman@wisc.edu}; Brian Yandell }

\keyword{hplot}
