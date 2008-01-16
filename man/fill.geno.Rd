\name{fill.geno}
\alias{fill.geno}

\title{Fill holes in genotype data}

\description{
  Replace the genotype data for a cross with a version imputed either
  by simulation with \code{\link[qtl]{sim.geno}} or by the Viterbi
  algorithm with \code{\link[qtl]{argmax.geno}}.
}

\usage{
fill.geno(cross, method=c("imp","argmax"), error.prob=0.0001,
          map.function=c("haldane","kosambi","c-f","morgan"))
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{method}{Indicates whether to impute using a single simulation
   replicate from \code{\link[qtl]{sim.geno}} or using the Viterbi
   algorithm, as implemented in \code{\link[qtl]{argmax.geno}}.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
 of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi or
   Carter-Falconer map function when converting genetic distances into
   recombination fractions.} 
}

\details{
  This function is written so that one may perform rough genome scans by
  marker regression without having to drop individuals with missing
  genotype data.  \bold{We must caution the user that little trust
  should be placed in the results}. 
}

\value{
  The input \code{cross} object with the genotype data replaced by an
  imputed version.  Any intermediate calculations (such as is produced
  by \code{\link[qtl]{calc.genoprob}}, \code{\link[qtl]{argmax.geno}}
  and \code{\link[qtl]{sim.geno}}) are removed. 
}

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper,chr=18:19)}
out.mr <- scantwo(fill.geno(hyper,method="argmax"), method="mr")
plot(out.mr)
}

\seealso{ \code{\link[qtl]{sim.geno}},
          \code{\link[qtl]{argmax.geno}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{utilities}
