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
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)

output <- comparegeno(listeria)

# image of the proportions
n.ind <- nind(listeria)
image(1:n.ind, 1:n.ind, output, col=gray((0:99)/99),
      breaks=seq(0,1,len=101))

# histogram
hist(output, breaks=150, prob=TRUE,
     xlab="Proportion of matching genotypes")
rug(output)
}

\seealso{ \code{\link{nmissing}} }

\keyword{utilities}
