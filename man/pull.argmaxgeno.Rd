\name{pull.argmaxgeno}
\alias{pull.argmaxgeno}

\title{Pull out the results of the Viterbi algorithm from a cross}

\description{
  Pull out the results of \code{\link{argmax.geno}} from a cross as a matrix.
}

\usage{pull.argmaxgeno(cross, chr, include.pos.info=FALSE, rotate=FALSE)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{include.pos.info}{If TRUE, include columns with marker name,
 chromosmoe ID, and cM position.  (If \code{include.pos.info=TRUE}, we
 take \code{rotate=TRUE}.)}
 \item{rotate}{If TRUE, return matrix with individuals as columns and
 positions as rows.  If FALSE, rows correspond to individuals.}
}

\value{
  A matrix containing numeric indicators of the inferred genotypes.
  Multiple chromosomes are pasted together.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
listeria <- argmax.geno(listeria, step=1, stepwidth="max")
amg <- pull.argmaxgeno(listeria, chr=c(5,13), include.pos.info=TRUE, rotate=TRUE)
amg[1:5,1:10]
}

\seealso{ \code{\link{pull.geno}}, \code{\link{pull.genoprob}},
  \code{\link{pull.draws}}, \code{\link{argmax.geno}} }

\keyword{utilities}
