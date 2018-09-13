\name{pull.geno}
\alias{pull.geno}

\title{Pull out the genotype data from a cross}

\description{
  Pull out the genotype data from a cross object, as a single big matrix.
}

\usage{pull.geno(cross, chr)}

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
  A matrix of size n.ind x tot.mar.  The raw genotype data in the input
  cross object, with the chromosomes pasted together.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
dat <- pull.geno(listeria)

# image of the genotype data
image(1:ncol(dat),1:nrow(dat),t(dat),ylab="Individuals",xlab="Markers",
      col=c("red","yellow","blue","green","violet"))
abline(v=cumsum(c(0,nmar(listeria)))+0.5)
abline(h=nrow(dat)+0.5)
}

\seealso{ \code{\link{pull.pheno}}, \code{\link{pull.map}}
  \code{\link{pull.draws}}, \code{\link{pull.genoprob}},
  \code{\link{pull.argmaxgeno}} }

\keyword{utilities}
