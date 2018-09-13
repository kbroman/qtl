\name{condense.scantwo}
\alias{condense.scantwo}

\title{Condense the output from a 2-d genome scan}

\description{
  Produces a very condensed version of the output of \code{\link{scantwo}}.
}

\usage{
\method{condense}{scantwo}(object)
}

\arguments{
 \item{object}{An object of class \code{scantwo}, the output of the
   function \code{\link{scantwo}}.}
}

\value{
  An object of class \code{scantwocondensed}, containing just the
  maximum full, additive and interactive LOD scores, and the positions
  where they occured, on each pair of chromosomes.
}

\details{
  This produces a very reduced version of the output of
  \code{\link{scantwo}}, for which a summary may still be created
  via \code{\link{summary.scantwo}}, though plots can no longer be
  made. 
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"))}
fake.f2 <- calc.genoprob(fake.f2)

out2 <- scantwo(fake.f2, method="hk")

out2c <- condense(out2)
summary(out2c, allpairs=FALSE)
max(out2c)
}

\seealso{ \code{\link{scantwo}}, \code{\link{summary.scantwo}},
  \code{\link{max.scantwo}} }

\keyword{print}
