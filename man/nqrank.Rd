\name{nqrank}
\alias{nqrank}

\title{Transform a vector of quantitative values to the corresponding
  normal quantiles}

\description{
  Transform a vector of quantitative values to the corresponding normal
  quantiles (preserving the mean and SD).
}

\usage{
nqrank(x, jitter)
}

\arguments{
  \item{x}{A numeric vector}
  \item{jitter}{If TRUE, randomly jitter the values to break ties.}
}

\value{
  A numeric vector; the input \code{x} is converted to ranks and then to
  normal quantiles.
}

\author{
Karl W Broman, \email{broman@wisc.edu}
}

\seealso{ \code{\link{rank}}, \code{\link{qnorm}}, \code{transformPheno} }

\examples{
data(hyper)

hyper <- transformPheno(hyper, pheno.col=1, transf=nqrank)
}

\keyword{ utilities }
