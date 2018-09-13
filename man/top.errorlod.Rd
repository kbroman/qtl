\name{top.errorlod}
\alias{top.errorlod}

\title{List genotypes with large error LOD scores}

\description{
  Prints those genotypes with error LOD scores above a specified
  cutoff.
}

\usage{
top.errorlod(cross, chr, cutoff=4, msg=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{cutoff}{Only those genotypes with error LOD scores above this
   cutoff will be listed.}
 \item{msg}{If TRUE, print a message if there are no apparent errors.}
}

\value{
  A data.frame with 4 columns, whose rows correspond to the genotypes
  that are possibly in error.  The four columns give the chromosome
  number, individual number, marker name, and error LOD score.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,19))}
# Calculate error LOD scores
hyper <- calc.errorlod(hyper,error.prob=0.01)

# Print those above a specified cutoff
top.errorlod(hyper,cutoff=4)
}

\seealso{ \code{\link{calc.errorlod}}, \code{\link{plotGeno}}, \code{\link{plotErrorlod}} }


\keyword{print}
