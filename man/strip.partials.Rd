\name{strip.partials}
\alias{strip.partials}

\title{Strip partially informative genotypes}

\description{
  Replace all partially informative genotypes (e.g., dominant markers in
  an intercross) with missing values.
}

\usage{strip.partials(cross, verbose=TRUE)}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{verbose}{If TRUE, print the number of genotypes removed.}
}

\value{
  The same class \code{cross} object as in the input, but with partially
  informative genotypes made missing.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
sum(nmissing(listeria))
listeria <- strip.partials(listeria)
sum(nmissing(listeria))
}

\seealso{ \code{\link{plotMissing}}, \code{\link{plotInfo}} }

\keyword{manip}
