\name{plot.scanoneboot}
\alias{plot.scanoneboot}

\title{Plot results of bootstrap for QTL position}

\description{
  Plot a histogram of the results of a nonparametric bootstrap to assess
  uncertainty in QTL position.
}

\usage{
\method{plot}{scanoneboot}(x, \dots)
}			    

\arguments{
  \item{x}{An object of class \code{"scanoneboot"}, as output by
    \code{\link{scanoneboot}}.}
  \item{\dots}{Passed to the function \code{\link[graphics]{hist}} when it
    is called.}
}
\value{None.}

\details{
  The function plots a histogram of the bootstrap results obtained by
  \code{\link{scanoneboot}}.  Genetic marker locations are
  displayed by vertical lines at the bottom of the plot.
}

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2, step=1)

\dontrun{out.boot <- scanoneboot(fake.f2, chr=13, method="hk")}
\dontshow{out.boot <- scanoneboot(fake.f2, chr=13, method="hk", n.boot=50)}

summary(out.boot)
plot(out.boot)
}

\seealso{ \code{\link{scanone}},
  \code{\link{summary.scanoneboot}} } 

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
