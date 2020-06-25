\name{effectscan}
\alias{effectscan}

\title{Plot estimated QTL effects across the whole genome}

\description{
  This function is used to plot the estimated QTL effects along
  selected chromosomes. For a backcross, there will be only one line,
  representing the additive effect. For an intercross, there will be two
  lines, representing the additive and dominance effects.
}

\usage{
effectscan(cross, pheno.col=1, chr, get.se=FALSE, draw=TRUE,
           gap=25, ylim, mtick=c("line","triangle"),
           add.legend=TRUE, alternate.chrid=FALSE, \dots)
}

\arguments{
  \item{cross}{An object of class \code{cross}.}
  \item{pheno.col}{Column number in the phenotype matrix which to be
    drawn in the plot.  One may also give a character string matching
   a phenotype name.}
  \item{chr}{Optional vector indicating the chromosomes to be drawn in
  the plot. This should be a vector of character strings referring to
  chromosomes by name; numeric values are converted to strings.  Refer
  to chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
  \item{get.se}{If TRUE, estimated standard errors are calculated.}
  \item{draw}{If TRUE, draw the figure.}
  \item{gap}{Gap separating chromosomes (in cM).}
  \item{ylim}{Y-axis limits (optional).}
  \item{mtick}{Tick mark type for markers.}
  \item{add.legend}{If TRUE, add a legend.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
     axis labels, so that they may be more easily distinguished.}
  \item{\dots}{Passed to the function \code{\link[base]{plot}} when it
    is called.}
}
\value{
  The results are returned silently, as an object of class
  \code{"effectscan"}, which is the same as the form returned by the
  function \code{\link{scanone}}, though with estimated effects
  where LOD scores might be.  That is, it is a data frame with the first
  two columns being chromosome ID and position (in cM), and subsequent
  columns being estimated effects, and (if \code{get.se=TRUE}) standard
  errors.
}

\details{
  The results of \code{\link{sim.geno}} are required for taking
  account of missing genotype information.

  For a backcross, the additive effect is estimated as the difference
  between the phenotypic averages for heterozygotes and homozygotes.

  For recombinant inbred lines, the additive effect is estimated as half
  the difference between the phenotypic averages for the two homozygotes.

  For an intercross, the additive and dominance effects are estimated
  from linear regression on \eqn{a} and \eqn{d} with \eqn{a} = -1, 0, 1,
  for the AA, AB and BB genotypes, respectively, and \eqn{d} = 0, 1, 0,
  for the AA, AB and BB genotypes, respectively.

  As usual, the X chromosome is a bit more complicated.  We estimate
  separate additive effects for the two sexes, and for the two
  directions within females.

  There is an internal function \code{plot.effectscan} that creates
  the actual plot by calling \code{plot.scanone}.  In the case
  \code{get.se=TRUE}, colored regions indicate \eqn{\pm}{+/-} 1 SE.
}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(1,13,"X"))}
fake.f2 <- sim.geno(fake.f2, step=2.5, n.draws=16)

# allelic effect on whole genome
effectscan(fake.f2)

# on chromosome 13, include standard errors
effectscan(fake.f2, chr="13", mtick="triangle", get.se=TRUE)
}

\references{
  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Karl W. Broman, \email{broman@wisc.edu}}

\seealso{ \code{\link{effectplot}}, \code{\link{plotPXG}},
  \code{\link{sim.geno}} }

\keyword{hplot}
