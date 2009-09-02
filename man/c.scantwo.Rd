\name{c.scantwo}
\alias{c.scantwo}
\alias{cbind.scantwo}

\title{Combine columns from multiple scantwo results}

\description{
  Concatenate the columns from different runs of
  \code{\link[qtl]{scantwo}}.
}

\usage{
\method{c}{scantwo}(\dots)
\method{cbind}{scantwo}(\dots)
}

\arguments{
 \item{\dots}{A set of objects of class \code{scantwo}. These are
   the results from  \code{\link[qtl]{scantwo}} (with \code{n.perm=0}),
   generally run with different phenotypes or methods.  All must 
   conform with each other, meaning that
   \code{\link[qtl]{calc.genoprob}} and/or \code{\link[qtl]{sim.geno}}
   were run with the same values for \code{step} and \code{off.end} and
   with data having the same genetic map.}
}

\value{The concatenated input, as a \code{scantwo} object.}

\details{
  The aim of this function is to concatenate the results from multiple
  runs \code{\link[qtl]{scantwo}}, generally for
  different phenotypes and/or methods.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.bc)
fake.bc <- calc.genoprob(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=18:19)}
out2a <- scantwo(fake.bc, method="hk")
out2b <- scantwo(fake.bc, pheno.col=2, method="hk")

out2 <- c(out2a, out2b)
}

\seealso{ \code{\link[qtl]{summary.scantwo}},
  \code{\link[qtl]{scantwo}}, \code{\link[qtl]{c.scanone}} }

\keyword{manip}
