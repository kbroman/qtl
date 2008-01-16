\name{geno.image}
\alias{geno.image}

\title{Plot grid of genotype data}

\description{
  Plot a grid showing which the genotype data in a cross.
}

\usage{
geno.image(x, chr, reorder=FALSE, main="Genotype data",
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
  Uses \code{\link[graphics]{image}} to plot a grid with the genotype
  data.  The genotypes AA, AB, BB are displayed in the colors red, blue,
  and green, respectively.  In an intercross, if there are genotypes
  ``not BB'' and ``not AA'', these are displayed in purple and orange,
  respectively.  White pixels indicate missing data.
}

\value{None.}

\examples{
data(listeria)
geno.image(listeria)
}

\seealso{ \code{\link[qtl]{plot.cross}},
  \code{\link[qtl]{plot.missing}}, \code{\link[qtl]{plot.geno}},
  \code{\link[graphics]{image}} }


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
