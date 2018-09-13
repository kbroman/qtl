\name{calc.errorlod}
\alias{calc.errorlod}

\title{Identify likely genotyping errors}

\description{
  Calculates a LOD score for each genotype, measuring the evidence for
  genotyping errors.
}

\usage{
calc.errorlod(cross, error.prob=0.01,
              map.function=c("haldane","kosambi","c-f","morgan"),
              version=c("new","old"))
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype)}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
 \item{version}{Specifies whether to use the original version of this
   function or the current (preferred) version.}
}

\details{
  Calculates, for each individual at each marker, a LOD score
  measuring the strength of evidence for a genotyping error, as
  described by Lincoln and Lander (1992).

  In the latest version, evidence for a genotype being in
  error is considered assuming that all other genotypes (for that
  individual, on that chromosome) are correct.  The argument
  \code{version} allows one to specify whether this new version is used,
  or whether the original (old) version of the calculation is
  performed.

  Note that values below 4 are generally not interesting.  Also note
  that if markers are extremely tightly linked, \emph{recombination
  events} can give large error LOD scores.  The error LOD scores should
  not be trusted blindly, but should be viewed as a tool for identifying
  genotypes deserving further study.

  Use \code{\link{top.errorlod}} to print all genotypes with error
  LOD scores above a specified threshold,
  \code{\link{plotErrorlod}} to plot the error LOD scores for
  specified chromosomes, and \code{\link{plotGeno}} to view the
  observed genotype data with likely errors flagged.
}

\value{
  The input \code{cross} object is returned with a component,
  \code{errorlod}, added to each component of \code{cross$geno}.  The
  \code{errorlod} component is a matrix of size (n.ind x n.mar).  An
  attribute \code{"error.prob"} is set to the value of the corresponding
  argument, for later reference.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\references{
  Lincoln, S. E. and Lander, E. S. (1992) Systematic detection of
  errors in genetic linkage data.  \emph{Genomics} \bold{14}, 604--610.
}

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper,chr=c(1,19))}
hyper <- calc.errorlod(hyper,error.prob=0.01)

# print those above a specified cutoff
top.errorlod(hyper, cutoff=4)

# plot genotype data, flagging genotypes with error LOD > cutoff
plotGeno(hyper, chr=1, ind=160:200, cutoff=7, min.sep=2)
}

\seealso{ \code{\link{plotErrorlod}},
  \code{\link{top.errorlod}}, \code{\link{cleanGeno}} }


\keyword{utilities}
