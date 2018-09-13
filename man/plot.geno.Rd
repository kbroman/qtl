\name{plotGeno}
\alias{plotGeno}

\title{Plot observed genotypes, flagging likely errors}

\description{
  Plot the genotypes on a particular chromosome for a set of
  individuals, flagging likely errors.
}

\usage{
plotGeno(x, chr, ind, include.xo=TRUE, horizontal=TRUE,
          cutoff=4, min.sep=2, cex=1.2, \dots)
}

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link{read.cross}} for details. }
  \item{chr}{The chromosome to plot.  Only one chromosome is
  allowed. (This should be a character string referring to the
  chromosomes by name.)}
  \item{ind}{Vector of individuals to plot (passed to \code{\link{subset.cross}}).  If missing, all individuals
    are plotted.}
  \item{include.xo}{If TRUE, plot X's in intervals having a
    crossover. Not available for a 4-way cross.}
  \item{horizontal}{If TRUE, chromosomes are plotted horizontally.}
  \item{cutoff}{Genotypes with error LOD scores above this value are
    flagged as possible errors.}
  \item{min.sep}{Markers separated by less than this value (as a
    percent of the chromosome length) are pulled apart, so that they may
    be distinguished in the picture.}
  \item{cex}{Character expansion for the size of points in the plot.
    Larger numbers give larger points; see \code{\link[graphics]{par}}.}
  \item{\dots}{Ignored at this point.}
}

\details{
  Plots the genotypes for a set of individuals.  Likely errors are
  indicated by red squares.  In a backcross, genotypes AA and AB are
  indicated by white and black circles, respectively.  In an intercross,
  genotypes AA, AB and BB are indicated by white, gray, and black
  circles, respectively, and the partially missing genotypes "not BB"
  (D in mapmaker) and "not AA" (C in mapmaker) are indicated by green
  and orange circles, respectively.

  For the X chromosome in a backcross or intercross, hemizygous males
  are plotted as if they were homozygous (that is, with white and black
  circles).

  For a 4-way cross, two lines are plotted for each individual.  The
  left or upper line indicates the allele A (white) or B (black); the
  right or lower line indicates the allele C (white) or D (black).  For
  the case that genotype is known to be only AC/BD or AD/BC, we use
  green and orange, respectively.

}

\value{None.}

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,19))}
# Calculate error LOD scores
hyper <- calc.errorlod(hyper,error.prob=0.01)

# print those above a specified cutoff
top.errorlod(hyper,cutoff=4)

# plot genotype data, flagging genotypes with error LOD > cutoff
plotGeno(hyper, chr=1, ind=160:200, cutoff=7, min.sep=2)
}

\seealso{ \code{\link{calc.errorlod}},
  \code{\link{top.errorlod}}, \code{\link{subset.cross}} }


\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
