\name{fill.geno}
\alias{fill.geno}

\title{Fill holes in genotype data}

\description{
  Replace the genotype data for a cross with a version imputed either
  by simulation with \code{\link{sim.geno}}, by the Viterbi
  algorithm with \code{\link{argmax.geno}}, or simply filling in
  genotypes between markers that have matching genotypes.
}

\usage{
fill.geno(cross, method=c("imp","argmax", "no_dbl_XO", "maxmarginal"),
          error.prob=0.0001,
          map.function=c("haldane","kosambi","c-f","morgan"),
          min.prob=0.95)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{method}{Indicates whether to impute using a single simulation
   replicate from \code{\link{sim.geno}}, using the Viterbi
   algorithm, as implemented in \code{\link{argmax.geno}}, by simply
   filling in missing genotypes between markers with matching genotypes,
   or by choosing (at each marker) the genotype with maximal marginal probability.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
 of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or
   Carter-Falconer map function when converting genetic distances into
   recombination fractions.} 
 \item{min.prob}{For \code{method="maxmarginal"}, genotypes with
 probability greater than this value will be imputed; those less than
 this value will be made missing.}
}

\details{
  This function is written so that one may perform rough genome scans by
  marker regression without having to drop individuals with missing
  genotype data.  \bold{We must caution the user that little trust
    should be placed in the results}.

  With \code{method="imp"}, a single random imputation is performed,
  using \code{\link{sim.geno}}.

  With \code{method="argmax"}, for each individual the most probable
  sequence of genotypes, given the observed data (via
  \code{argmax.geno}), is used.

  With \code{method="no_dbl_XO"}, non-recombinant intervals are filled
  in; recombinant intervals are left missing. For example, a sequence of
  genotypes like \code{A---A---H---H---A} (with \code{A} and \code{H}
  corresponding to genotypes AA and AB, respectively, and with \code{-}
  being a missing value) will be filled in as
  \code{AAAAA---HHHHH---A}.

  With \code{method="maxmarginal"}, the conditional genotype
  probabilities are calculated with \code{calc.genoprob}, and then at
  each marker, the most probable genotype is determined. This is taken
  as the imputed genotype if it has probability greater than
  \code{min.prob}; otherwise it is made missing.

  With \code{method="no_dbl_XO"} and \code{method="maxmarginal"}, 
  some missing genotypes likely remain. With
  \code{method="maxmarginal"}, some observed genotypes may be made
  missing.
}

\value{
  The input \code{cross} object with the genotype data replaced by an
  imputed version.  Any intermediate calculations (such as is produced
  by \code{\link{calc.genoprob}}, \code{\link{argmax.geno}}
  and \code{\link{sim.geno}}) are removed. 
}

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper,chr=18:19)}
out.mr <- scantwo(fill.geno(hyper,method="argmax"), method="mr")
plot(out.mr)
}

\seealso{ \code{\link{sim.geno}},
          \code{\link{argmax.geno}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{utilities}
