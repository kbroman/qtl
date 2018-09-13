\name{plotInfo}
\alias{plotInfo}

\title{Plot the proportion of missing genotype information}

\description{
  Plot a measure of the proportion of missing information in the
  genotype data.
}

\usage{
plotInfo(x, chr, method=c("entropy","variance","both"), step=1,
          off.end=0, error.prob=0.001,
          map.function=c("haldane","kosambi","c-f","morgan"),
          alternate.chrid=FALSE, fourwaycross=c("all", "AB", "CD"),
          include.genofreq=FALSE, \dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to plot.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{method}{Indicates whether to plot the entropy version of the
   information, the variance version, or both.}
 \item{step}{Maximum distance (in cM) between positions at which the
   missing information is calculated, though for \code{step=0},
   it is are calculated only at the marker locations.}
 \item{off.end}{Distance (in cM) past the terminal markers on each
   chromosome to which the genotype probability calculations will be
   carried.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or
   Carter-Falconer map function when converting genetic distances into
   recombination fractions.}
 \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
    axis labels, so that they may be more easily distinguished.}
  \item{fourwaycross}{For a phase-known four-way cross, measure missing
 genotype information overall (\code{"all"}), or just for the alleles
 from the first parent (\code{"AB"}) or from the second parent (\code{"CD"}).}
  \item{include.genofreq}{If TRUE, estimated genotype frequencies (from
 the results of
 \code{\link{calc.genoprob}} averaged across the individuals) are
 included as additional columns in the output.}
 \item{\dots}{Passed to \code{\link{plot.scanone}}.}
}

\details{
  The entropy version of the missing information: for a single
  individual at a single genomic position, we measure the missing
  information as \eqn{H = \sum_g p_g \log p_g / \log n}{H = sum p[g] log
  p[g] / log n}, where \eqn{p_g}{p[g]} is the probability of the
  genotype \eqn{g}, and \eqn{n} is the number of possible genotypes,
  defining \eqn{0 \log 0 = 0}{0 log 0 = 0}.  This takes values between 0
  and 1, assuming the value 1 when the genotypes (given the marker data)
  are equally likely and 0 when the genotypes are completely determined.
  We calculate the missing information at a particular position as the
  average of \eqn{H} across individuals.  For an intercross, we don't
  scale by \eqn{\log n} but by the entropy in the case of genotype
  probabilities (1/4, 1/2, 1/4).

  The variance version of the missing information: we calculate the
  average, across individuals, of the variance of the genotype
  distribution (conditional on the observed marker data) at a particular
  locus, and scale by the maximum such variance.

  Calculations are done in C (for the sake of speed in the presence of
  little thought about programming efficiency) and the plot is created
  by a call to \code{\link{plot.scanone}}.

  Note that \code{\link{summary.scanone}} may be used to display
  the maximum missing information on each chromosome.
}

\value{
  An object with class \code{scanone}: a data.frame with columns the
  chromosome IDs and cM positions followed by the entropy and/or
  variance version of the missing information.
}

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper,chr=1:4)}
plotInfo(hyper,chr=c(1,4))

# save the results and view maximum missing info on each chr
info <- plotInfo(hyper)
summary(info)

plotInfo(hyper, bandcol="gray70")
}

\seealso{ \code{\link{plot.scanone}},
  \code{\link{plotMissing}}, \code{\link{calc.genoprob}},
  \code{\link{geno.table}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
\keyword{univar}
