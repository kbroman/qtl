\name{plot.qtl}
\alias{plot.qtl}

\title{Plot QTL locations}

\description{
  Plot the locations of the QTL against a genetic map
}

\usage{
\method{plot}{qtl}(x, chr, horizontal=FALSE, shift=TRUE,
     show.marker.names=FALSE, alternate.chrid=FALSE, \dots)
}

\arguments{
  \item{x}{An object of class \code{"qtl"}, as produced by
    \code{\link[qtl]{makeqtl}}.} 
  \item{chr}{Vector specifying which chromosomes to plot. (The
    chromosomes must be specified by name.)  If missing, all chromosomes
    are plotted.}
  \item{horizontal}{Specifies whether the chromosomes should be
    plotted horizontally.}
  \item{shift}{If TRUE, shift the first marker on each chromosome to be
    at 0 cM.}
  \item{show.marker.names}{If TRUE, marker names are included.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome 
     axis labels, so that they may be more easily distinguished.}
  \item{\dots}{Passed to \code{\link[qtl]{plot.map}}.}
}

\value{None.}

\details{
  Creates a plot, via \code{\link[qtl]{plot.map}}, and indicates the
  locations of the QTL in the input QTL object, \code{x}.
}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c("1", "6", "13")
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- sim.geno(fake.f2, n.draws=2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="draws")
plot(qtl)
}

\seealso{ \code{\link[qtl]{plot.map}}, \code{\link[qtl]{makeqtl}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
