\name{pull.draws}
\alias{pull.draws}

\title{Pull out the genotype imputations from a cross}

\description{
  Pull out the results of \code{\link{sim.geno}} from a cross as an array.
}

\usage{pull.draws(cross, chr)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
}

\value{
  An array containing numeric indicators of the imputed genotypes.
  Multiple chromosomes are pasted together.  The dimensions are
  individuals by positions by imputations
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
listeria <- sim.geno(listeria, step=5, stepwidth="max", n.draws=8)
dr <- pull.draws(listeria, chr=c(5,13))
dr[1:20,1:10,1]
}

\seealso{ \code{\link{pull.geno}}, \code{\link{pull.genoprob}},
  \code{\link{pull.argmaxgeno}}, \code{\link{sim.geno}} }

\keyword{utilities}
