\name{sim.geno}
\alias{sim.geno}

\title{Simulate genotypes given observed marker data}

\description{
  Uses the hidden Markov model technology to simulate from the joint
  distribution Pr(g | O) where g is the underlying genotype vector and O
  is the observed multipoint marker data, with possible allowance for
  genotyping errors.
}

\usage{
sim.geno(cross, n.draws=16, step=0, off.end=0, error.prob=0.0001,
         map.function=c("haldane","kosambi","c-f","morgan"),
         stepwidth=c("fixed", "variable", "max"))
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{n.draws}{Number of simulation replicates to perform.}
 \item{step}{Maximum distance (in cM) between positions at which the
   simulated genotypes will be drawn, though for \code{step=0},
   genotypes are drawn only at the marker locations.}
 \item{off.end}{Distance (in cM) past the terminal markers on each
   chromosome to which the genotype simulations will be carried.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
 \item{stepwidth}{Indicates whether the intermediate points should with
   fixed or variable step sizes.  We recommend using
   \code{"fixed"}; \code{"variable"} was included for the qtlbim
   package (\url{https://cran.r-project.org/src/contrib/Archive/qtlbim/}).  The \code{"max"}
   option inserts the minimal number of intermediate points so that the
   maximum distance between points is \code{step}.}
}

\details{
  After performing the forward-backward equations, we draw from
  \eqn{Pr(g_1 = v | O)}{Pr(g[1] = v | O)} and then \eqn{Pr(g_{k+1} = v |
    O, g_k = u)}{Pr(g[k+1] = v | O, g[k] = u)}.

  In the case of the 4-way cross, with a sex-specific map, we assume a
  constant ratio of female:male recombination rates within the
  inter-marker intervals.
}

\value{
  The input \code{cross} object is returned with a component,
  \code{draws}, added to each component of \code{cross$geno}.
  This is an array of size [n.ind x n.pos x n.draws] where n.pos is
  the number of positions at which the simulations were performed and
  n.draws is the number of replicates.  Attributes \code{"error.prob"},
  \code{"step"}, and \code{"off.end"} are set to the values of the
  corresponding arguments, for later reference.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
fake.f2 <- sim.geno(fake.f2, step=2, n.draws=8)
}

\seealso{ \code{\link{calc.genoprob}},
  \code{\link{argmax.geno}} }

\keyword{utilities}
