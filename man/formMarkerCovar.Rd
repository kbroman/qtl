\name{formMarkerCovar}
\alias{formMarkerCovar}

\title{Create matrix of marker covariates for QTL analysis}

\description{
  Pull out a matrix of genotypes or genotype probabilities to use
  markers as covariates in QTL analysis.
}

\usage{formMarkerCovar(cross, markers, method=c("prob", "imp", "argmax"), \dots)}


\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{markers}{A vector of character strings of marker or pseudomarker
   names. Pseudomarker names may be of the form \code{"5@21.5"} (for chr 5
   at 21.5 cM), but then all names must be of this form.}
 \item{method}{If \code{method="prob"}, the genotype probabilities from
   \code{\link{calc.genoprob}} are used; otherwise we use
   \code{\link{fill.geno}} to impute missing data, with this method.}
 \item{\dots}{Passed to \code{\link{fill.geno}}, if necessary.}
}

\value{
  A matrix containing genotype probabilities or genotype indicators,
  suitable for use as covariates in \code{\link{scanone}}.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
hyper <- calc.genoprob(hyper, step=0)
peakMarker <- "D4Mit164"
X <- formMarkerCovar(hyper, peakMarker)
\dontshow{hyper <- hyper[c(1,4),]}
out <- scanone(hyper, addcovar=X)
}

\seealso{ \code{\link{pull.geno}}, \code{\link{pull.genoprob}},
  \code{\link{fill.geno}}, \code{\link{scanone}} }

\keyword{utilities}
