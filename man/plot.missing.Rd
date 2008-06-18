\name{plot.missing}
\alias{plot.missing}

\title{Plot grid of missing genotypes}

\description{
  Plot a grid showing which genotypes are missing.
}

\usage{
plot.missing(x, chr, reorder=FALSE, main="Missing genotypes",
             alternate.chrid=FALSE, \dots)
}			    

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link[qtl]{read.cross}} for details. }
  \item{chr}{Vector of chromosomes to plot.  If missing, plot all
    chromosomes.}
  \item{reorder}{Specify whether to reorder individuals according to
  their phenotypes.
  \tabular{ll}{
    FALSE \tab Don't reorder \cr
    TRUE  \tab Reorder according to the sum of the phenotypes \cr
    n \tab Reorder according to phenotype n}}
  \item{main}{Title to place on plot.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome 
     axis labels, so that they may be more easily distinguished.}
  \item{\dots}{Ignored at this point.}
}

\details{
  Uses \code{\link[graphics]{image}} to plot a grid with black pixels where the
  genotypes are missing.  For intercross and 4-way cross data, gray
  pixels are plotted for the partially missing genotypes (for example,
  "not AA").
}

\value{None.}

\examples{
data(fake.f2)
plot.missing(fake.f2)
}

\seealso{ \code{\link[qtl]{plot.cross}},  \code{\link[qtl]{geno.image}},
  \code{\link[graphics]{image}} }


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
