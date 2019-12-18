\name{comparegeno}
\alias{comparegeno}

\title{Compare individuals' genotype data}

\description{
  Count proportion of matching genotypes between all pairs of
  individuals, to look for unusually closely related individuals.
}

\usage{comparegeno(cross, what=c("proportion","number","both"))}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{what}{Indicates whether to return the proportion or number of matching
genotypes (or both).}
}

\value{
  A matrix whose (i,j)th element is the proportion or number of matching
  genotypes for individuals i and j.

  If called with \code{what="both"}, the lower triangle contains the
  proportion and the upper triangle contains the number.

  If called with \code{what="proportion"}, the diagonal contains missing
  values.  Otherwise, the diagonal contains the number of typed markers
  for each individual.

  The output is given class \code{"comparegeno"} so that appropriate
  \code{summary} and \code{plot} functions may be used.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)

cg <- comparegeno(listeria)

summary(cg, 0.7)
plot(cg)
}

\seealso{ \code{\link{nmissing}}, \code{\link{summary.comparegeno}},
          \code{\link{plot.comparegeno}} }

\keyword{utilities}
