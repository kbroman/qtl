\name{pull.genoprob}
\alias{pull.genoprob}

\title{Pull out the genotype probabilities from a cross}

\description{
  Pull out the results of \code{\link{calc.genoprob}} from a cross as a matrix.
}

\usage{pull.genoprob(cross, chr, omit.first.prob=FALSE,
              include.pos.info=FALSE, rotate=FALSE)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{omit.first.prob}{If TRUE, omit the probabilities for the first
 genotype at each position (since they sum to 1).}
 \item{include.pos.info}{If TRUE, include columns with marker name,
 genotype, chromosome ID, and cM position.  (If
 \code{include.pos.info=TRUE}, we take \code{rotate=TRUE}.)}
 \item{rotate}{If TRUE, return matrix with individuals as columns and
 positions/genotypes as rows.  If FALSE, rows correspond to individuals.}
}

\value{
  A matrix containing genotype probabilities.
  Multiple chromosomes and the multiple genotypes at each position are
  pasted together. 
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
listeria <- calc.genoprob(listeria, step=1, stepwidth="max")
pr <- pull.genoprob(listeria, chr=c(5,13), omit.first.prob=TRUE, include.pos.info=TRUE, rotate=TRUE)
pr[1:5,1:10]
}

\seealso{ \code{\link{pull.geno}}, \code{\link{pull.argmaxgeno}},
  \code{\link{pull.draws}}, \code{\link{calc.genoprob}} }

\keyword{utilities}
