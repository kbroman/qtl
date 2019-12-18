\name{summary.comparegeno}
\alias{summary.comparegeno}

\title{Print pairs of individuals with similar genotype data.}

\description{
  Prints a summary the output from \code{\link{comparegeno}} that
  includes pairs of individuals whose proportion of matching genotypes
  is above a chosen threshold.
}

\usage{
\method{summary}{comparegeno}(object, thresh=0.9, \dots)
}

\arguments{
 \item{object}{An object of class \code{comparegeno}, the output of the
   function \code{\link{comparegeno}}.}
 \item{thresh}{Threshold on the proportion of matching genotypes.}
 \item{\dots}{Ignored at this point.}
}

\value{
  A data frame with each row being a pair of individuals and columns
  including the individual identifiers (via \code{getid}, or just as
  numeric indexes) along with the proportion of matching genotypes.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
cg <- comparegeno(fake.f2)
summary(cg, 0.7)
}

\seealso{ \code{\link{comparegeno}}, \code{\link{plot.comparegeno}} }

\keyword{print}
