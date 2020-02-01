\name{argmax.geno}
\alias{argmax.geno}

\title{Reconstruct underlying genotypes}

\description{
  Uses the Viterbi algorithm to identify the most likely sequence
  of underlying genotypes, given the observed multipoint marker data,
  with possible allowance for genotyping errors.
}

\usage{
argmax.geno(cross, step=0, off.end=0, error.prob=0.0001,
            map.function=c("haldane","kosambi","c-f","morgan"),
            stepwidth=c("fixed", "variable", "max"))
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{step}{Maximum distance (in cM) between positions at which the
   genotypes are reconstructed, though for \code{step=0}, genotypes
   are reconstructed only at the marker locations.}
 \item{off.end}{Distance (in cM) past the terminal markers on each
   chromosome to which the genotype reconstructions will be carried.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer or Morgan map function when converting genetic
   distances into recombination fractions.}
 \item{stepwidth}{Indicates whether the intermediate points should with
   fixed or variable step sizes.  We recommend using
   \code{"fixed"}; \code{"variable"} was included for the qtlbim
   package (\url{https://cran.r-project.org/src/contrib/Archive/qtlbim}).  The \code{"max"}
   option inserts the minimal number of intermediate points so that the
   maximum distance between points is \code{step}.}
}

\details{
  We use the Viterbi algorithm to calculate
  \eqn{\arg \max_v \Pr(g = v | O)}{arg max_v Pr(g = v | O)} where
  \eqn{g} is the underlying sequence of genotypes and \eqn{O} is the
  observed marker genotypes.

  This is done by calculating
  \eqn{\gamma_k(v_k) = \max_{v_1, \ldots, v_{k-1}} \Pr(g_1 = v_1,
    \ldots, g_k = v_k, O_1, \ldots, O_k)}{%
    Q[k](v[k]) = max{v[1], \ldots, v[k-1]} Pr(g[1] = v[1],
    \ldots, g[k] = v[k], O[1], \ldots, O[k])}
  for \eqn{k = 1, \ldots, n} and then tracing back through the
  sequence.
}

\value{
  The input \code{cross} object is returned with a component,
  \code{argmax}, added to each component of \code{cross$geno}.
  The \code{argmax} component is a matrix of size [n.ind x n.pos], where
  n.pos is the
  number of positions at which the reconstructed genotypes were obtained,
  containing the most likely sequences of underlying genotypes.
  Attributes \code{"error.prob"}, \code{"step"}, and \code{"off.end"}
  are set to the values of the corresponding arguments, for later
  reference.
}

\section{Warning}{
  The Viterbi algorithm can behave badly when \code{step} is small but
  positive. One may observe quite different results for different values
  of \code{step}.

  The problem is that, in the presence of data like \code{A----H},  the
  sequences \code{AAAAAA} and \code{HHHHHH} may be more likely than any
  one of the sequences \code{AAAAAH}, \code{AAAAHH}, \code{AAAHHH},
  \code{AAHHHH}, \code{AHHHHH}, \code{AAAAAH}.  The Viterbi algorithm
  produces a single "most likely" sequence of underlying genotypes.
}


\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2,chr=18:19)
}fake.f2 <- argmax.geno(fake.f2, step=2, off.end=5, err=0.01)
}

\references{
  Lange, K. (1999) \emph{Numerical analysis for statisticians}.
  Springer-Verlag. Sec 23.3.

  Rabiner, L. R. (1989) A tutorial on hidden Markov models and selected
  applications in speech recognition.  \emph{Proceedings of the IEEE}
  \bold{77}, 257--286.
}


\seealso{ \code{\link{sim.geno}}, \code{\link{calc.genoprob}},
\code{\link{fill.geno}} }

\keyword{utilities}
