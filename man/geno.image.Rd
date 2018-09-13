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
    \code{\link{read.cross}} for details. }
  \item{chr}{Optional vector indicating the chromosomes to be drawn in
  the plot. This should be a vector of character strings referring to
  chromosomes by name; numeric values are converted to strings.  Refer
  to chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
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
  \item{\dots}{Passed to \code{\link[graphics]{image}}.}
}

\details{
  Uses \code{\link[graphics]{image}} to plot a grid with the genotype
  data.  The genotypes AA, AB, BB are displayed in the colors red, blue,
  and green, respectively.  In an intercross, if there are genotypes
  "not BB" and "not AA", these are displayed in purple and orange,
  respectively.  White pixels indicate missing data.
}

\value{None.}

\examples{
data(listeria)
geno.image(listeria)
}

\seealso{ \code{\link{plot.cross}},
  \code{\link{plotMissing}}, \code{\link{plotGeno}},
  \code{\link[graphics]{image}} }


\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
