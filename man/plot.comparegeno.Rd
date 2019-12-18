\name{plot.comparegeno}
\alias{plot.comparegeno}

\title{Plot genotype comparison}

\description{
  Plot the results of the comparison of all pairs of individuals'
  genotypes. A histogram of the proportion of matching genotypes, with
  tick marks at individual values below, via \code{\link[graphics]{rug}}.
}

\usage{
\method{plot}{comparegeno}(x, breaks=NULL, main="",
     xlab="Proportion matching genotypes", \dots)
}

\arguments{
  \item{x}{An object of class \code{"comparegeno"}, as produced by
    \code{\link{comparegeno}}.}
  \item{breaks}{Passed to \code{\link[graphics]{hist}}, with the default
    \eqn{2 \sqrt{n}} where \eqn{n} is the number of pairs of individuals.}
  \item{main}{Title for the plot.}
  \item{xlab}{x-axis label for the plot.}
  \item{\dots}{Passed to \code{\link{hist}}.}
}

\value{None.}

\details{
  Creates a histogram with \code{\link[graphics]{hist}} with ticks at
  individual values using \code{\link[graphics]{rug}}.
}

\examples{
data(fake.f2)
cg <- comparegeno(fake.f2)
plot(cg)
}

\seealso{ \code{\link{comparegeno}}, \code{\link{summary.comparegeno}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
