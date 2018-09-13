\name{calc.genoprob}
\alias{calc.genoprob}

\title{Calculate conditional genotype probabilities}

\description{
  Uses the hidden Markov model technology to calculate the
  probabilities of the true underlying genotypes given the observed
  multipoint marker data, with possible allowance for genotyping
  errors.
}

\usage{
calc.genoprob(cross, step=0, off.end=0, error.prob=0.0001,
              map.function=c("haldane","kosambi","c-f","morgan"),
              stepwidth=c("fixed", "variable", "max"))
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{step}{Maximum distance (in cM) between positions at which the
   genotype probabilities are calculated, though for \code{step = 0},
   probabilities are calculated only at the marker locations.}
 \item{off.end}{Distance (in cM) past the terminal markers on each
   chromosome to which the genotype probability calculations will be
   carried.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or
   Carter-Falconer map function when converting genetic distances into
   recombination fractions.}
 \item{stepwidth}{Indicates whether the intermediate points should with
   fixed or variable step sizes.  We recommend using
   \code{"fixed"}; \code{"variable"} is included for the qtlbim
   package (\url{http://www.ssg.uab.edu/qtlbim}).  The \code{"max"}
   option inserts the minimal number of intermediate points so that the
   maximum distance between points is \code{step}.}
}

\details{
  Let \eqn{O_k}{O[k]} denote the observed marker genotype at position
  \eqn{k}, and \eqn{g_k}{g[k]} denote the corresponding true underlying
  genotype.

  We use the forward-backward equations to calculate
  \eqn{\alpha_{kv} = \log Pr(O_1, \ldots, O_k, g_k = v)}{%
    a[k][v] = log Pr(O[1], \ldots, O[k], g[k] = v)}
  and
  \eqn{\beta_{kv} = \log Pr(O_{k+1}, \ldots, O_n | g_k = v)}{%
    b[k][v] = log Pr(O[k+1], \ldots, O[n] | g[k] = v)}

  We then obtain
  \eqn{Pr(g_k | O_1, \ldots, O_n) = \exp(\alpha_{kv} + \beta_{kv}) / s}{%
    Pr(g[k] | O[1], \ldots, O[n] = exp(a[k][v] + b[k][v]) / s}
  where
  \eqn{s = \sum_v \exp(\alpha_{kv} + \beta_{kv})}{%
    s = sum_v exp(a[k][v] + b[k][v])}

  In the case of the 4-way cross, with a sex-specific map, we assume a
  constant ratio of female:male recombination rates within the
  inter-marker intervals.
}

\value{
  The input \code{cross} object is returned with a component,
  \code{prob}, added to each component of \code{cross$geno}.
  \code{prob} is an array of size [n.ind x n.pos x n.gen] where n.pos is
  the number of positions at which the probabilities were calculated and
  n.gen = 3 for an intercross, = 2 for a backcross, and = 4 for a 4-way
  cross.  Attributes \code{"error.prob"}, \code{"step"},
  \code{"off.end"}, and \code{"map.function"} are set to the values of
  the corresponding arguments, for later reference (especially by the
  function \code{\link{calc.errorlod}}).
}

\references{
  Lange, K. (1999) \emph{Numerical analysis for statisticians}.
  Springer-Verlag. Sec 23.3.

  Rabiner, L. R. (1989) A tutorial on hidden Markov models and selected
  applications in speech recognition.  \emph{Proceedings of the IEEE}
  \bold{77}, 257--286.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2,chr=18:19)
}fake.f2 <- calc.genoprob(fake.f2, step=2, off.end=5)

data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc,chr=18:19)
}fake.bc <- calc.genoprob(fake.bc, step=0, off.end=0, err=0.01)
}

\seealso{ \code{\link{sim.geno}}, \code{\link{argmax.geno}},
         \code{\link{calc.errorlod}} }

\keyword{utilities}
