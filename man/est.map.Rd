\name{est.map}
\alias{est.map}

\title{Estimate genetic maps}

\description{
  Uses the Lander-Green algorithm (i.e., the hidden Markov model
  technology) to re-estimate the genetic map for an experimental cross.
}

\usage{
est.map(cross, chr, error.prob=0.0001,
        map.function=c("haldane","kosambi","c-f","morgan"),
        m=0, p=0, maxit=10000, tol=1e-6, sex.sp=TRUE,
        verbose=FALSE, omit.noninformative=TRUE, offset, n.cluster=1)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
   This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions. (Ignored if m > 0.)}
 \item{m}{Interference parameter for the chi-square model for
   interference; a non-negative integer, with m=0 corresponding to no
   interference. This may be used only for a backcross or intercross.}
 \item{p}{Proportion of chiasmata from the NI mechanism, in the Stahl
   model; p=0 gives a pure chi-square model.  This may be used only for
   a backcross or intercross.} 
 \item{maxit}{Maximum number of EM iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps; this is 
 used only for the 4-way cross.}
 \item{verbose}{If TRUE, print tracing information.}
 \item{omit.noninformative}{If TRUE, on each chromosome, omit individuals
   with fewer than two typed markers, since they are not informative for
   linkage.}
 \item{offset}{Defines the starting position
  for each chromosome.  If missing, we use the starting positions that
  are currently present in the input cross object.  This should be a
  single value (to be used for all chromosomes) or a vector with length
  equal to the number of chromosomes, defining individual starting
  positions for each chromosome.  For a sex-specific map (as in a 4-way
  cross), we use the same offset for both the male and female maps.}
 \item{n.cluster}{If the package \code{snow} is available 
   calculations for multiple chromosomes are run in parallel using this
   number of nodes.}
}

\details{
  By default, the map is estimated assuming no crossover interference,
  but a map function is used to derive the genetic distances (though, by
  default, the Haldane map function is used).

  For a backcross or intercross, inter-marker distances may be estimated
  using the Stahl model for crossover interference, of which the
  chi-square model is a special case.

  In the chi-square model, points are tossed down onto the four-strand
  bundle according to a Poisson process, and every \eqn{(m+1)}st point is a
  chiasma.  With the assumption of no chromatid interference, crossover
  locations on a random meiotic product are obtained by thinning the
  chiasma process.  The parameter \eqn{m} (a non-negative integer)
  governs the strength of crossover interference, with \eqn{m=0}
  corresponding to no interference.

  In the Stahl model, chiasmata on the four-strand bundle are a
  superposition of chiasmata from two mechanisms, one following a
  chi-square model and one exhibiting no interference.  An additional
  parameter, \eqn{p}, gives the proportion of chiasmata from the no
  interference mechanism.
}

\value{
  A \code{map} object; a list whose components (corresponding to
  chromosomes) are either vectors of marker positions (in cM) or
  matrices with two rows of sex-specific marker positions.
  The maximized log likelihood for each chromosome is saved as an
  attribute named \code{loglik}.  In the case that estimation was under
  an interference model (with m > 0), allowed only for a backcross, m
  and p are also included as attributes.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2,chr=18:19)}
newmap <- est.map(fake.f2)
logliks <- sapply(newmap, attr, "loglik")
plotMap(fake.f2, newmap)
fake.f2 <- replace.map(fake.f2, newmap)
}

\references{
  Armstrong, N. J., McPeek, M. J. and Speed, T. P. (2006) Incorporating
  interference into linkage analysis for experimental crosses.
  \emph{Biostatistics} \bold{7}, 374--386.
  
  Lander, E. S. and Green, P. (1987) Construction of multilocus genetic linkage
  maps in humans.  \emph{Proc. Natl. Acad. Sci. USA} \bold{84}, 2363--2367.

  Lange, K. (1999) \emph{Numerical analysis for statisticians}.
  Springer-Verlag. Sec 23.3.

  Rabiner, L. R. (1989) A tutorial on hidden Markov models and selected
  applications in speech recognition.  \emph{Proceedings of the IEEE}
  \bold{77}, 257--286.

  Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis of
  crossover interference using the chi-square model.  \emph{Genetics}
  \bold{139}, 1045--1056.
}  

\seealso{ \code{\link{plotMap}}, \code{\link{replace.map}},
  \code{\link{est.rf}}, \code{\link{fitstahl}} }

\keyword{utilities}
