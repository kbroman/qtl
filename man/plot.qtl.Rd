\name{plot.qtl}
\alias{plot.qtl}

\title{Plot QTL locations}

\description{
  Plot the locations of the QTL against a genetic map
}

\usage{
\method{plot}{qtl}(x, chr, horizontal=FALSE, shift=TRUE,
     show.marker.names=FALSE, alternate.chrid=FALSE, justdots=FALSE,
     col="red", \dots)
}

\arguments{
  \item{x}{An object of class \code{"qtl"}, as produced by
    \code{\link{makeqtl}}.} 
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
  \item{justdots}{If FALSE, just plot dots at the QTL, rather than
    arrows and QTL names.}
  \item{col}{Color used to plot indications of QTL}
  \item{\dots}{Passed to \code{\link{plotMap}}.}
}

\value{None.}

\details{
  Creates a plot, via \code{\link{plotMap}}, and indicates the
  locations of the QTL in the input QTL object, \code{x}.
}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c("1", "6", "13")
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")
plot(qtl)
plot(qtl, justdots=TRUE, col="seagreen")
}

\seealso{ \code{\link{plotMap}}, \code{\link{makeqtl}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
