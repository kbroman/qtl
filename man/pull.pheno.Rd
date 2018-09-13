\name{pull.pheno}
\alias{pull.pheno}

\title{Pull out phenotype data from a cross}

\description{
  Pull out selected phenotype data from a cross object, as a data frame
  or vector.
}

\usage{pull.pheno(cross, pheno.col)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{A vector specifying which phenotypes to keep or discard.
   This may be a logical vector, a numeric vector, or a vector
   of character strings (for the phenotype names). If missing, the
   entire set of phenotypes is output.} 
}

\value{
  A data.frame with columns specifying phenotypes and rows specifying individuals.
  If there is just one phenotype, a vector (rather than a data.frame) is returned.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
pull.pheno(listeria, "sex")
}

\seealso{ \code{\link{pull.geno}}, \code{\link{pull.map}} }

\keyword{utilities}
