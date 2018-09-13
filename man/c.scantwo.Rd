\name{c.scantwo}
\alias{c.scantwo}
\alias{cbind.scantwo}

\title{Combine columns from multiple scantwo results}

\description{
  Concatenate the columns from different runs of
  \code{\link{scantwo}}.
}

\usage{
\method{c}{scantwo}(\dots)
\method{cbind}{scantwo}(\dots)
}

\arguments{
 \item{\dots}{A set of objects of class \code{scantwo}. (This can also
 be a list of \code{scantwo} objects.) These are
   the results from  \code{\link{scantwo}} (with \code{n.perm=0}),
   generally run with different phenotypes or methods.  All must
   conform with each other, meaning that
   \code{\link{calc.genoprob}} and/or \code{\link{sim.geno}}
   were run with the same values for \code{step} and \code{off.end} and
   with data having the same genetic map.}
}

\value{The concatenated input, as a \code{scantwo} object.}

\details{
  The aim of this function is to concatenate the results from multiple
  runs \code{\link{scantwo}}, generally for
  different phenotypes and/or methods.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=18:19)}
out2a <- scantwo(fake.bc, method="hk")
out2b <- scantwo(fake.bc, pheno.col=2, method="hk")

out2 <- c(out2a, out2b)
}

\seealso{ \code{\link{summary.scantwo}},
  \code{\link{scantwo}}, \code{\link{c.scanone}} }

\keyword{manip}
