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
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector specifying for which chromosomes the error LOD
   scores should be inspected.}
 \item{cutoff}{Only those genotypes with error LOD scores above this
   cutoff will be listed.}
 \item{msg}{If TRUE, print a message if there are no apparent errors.}
}

\value{
  A data.frame with 4 columns, whose rows correspond to the genotypes
  that are possibly in error.  The four columns give the chromosome
  number, individual number, marker name, and error LOD score.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,19))}
# Calculate error LOD scores
hyper <- calc.errorlod(hyper,error.prob=0.01)

# Print those above a specified cutoff 
top.errorlod(hyper,cutoff=4)
}

\seealso{ \code{\link[qtl]{calc.errorlod}}, \code{\link[qtl]{plot.geno}} }


\keyword{print}
