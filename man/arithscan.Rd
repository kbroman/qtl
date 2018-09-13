\name{arithscan}
\alias{+.scanone}
\alias{-.scanone}
\alias{+.scantwo}
\alias{-.scantwo}

\title{Arithmetic operators for scanone and scantwo results}

\description{
  Add or subtract LOD scores in results
  from \code{\link{scanone}} or \code{\link{scantwo}}.
}

\usage{
\special{scan1+scan2}
\special{scan1-scan2}
}

\arguments{
  \item{scan1, scan2}{Genome scan results on the same set of chromosomes
    and markers, as output by \code{\link{scanone}} or \code{\link{scantwo}}.}
}

\details{
  This is used to calculate the sum or difference of LOD scores of two
  genome scan results.  It is particularly useful for calculating the LOD
  scores for QTL-by-covariate interactions (see the example, below).
  Note that the degrees of freedom are also added or subtracted.
}

\value{
  The same type of data structure as the input objects, with LOD scores
  added or subtracted.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=1:2)}
fake.bc <- calc.genoprob(fake.bc, step=2.5)

# covariates
ac <- pull.pheno(fake.bc, c("sex","age"))
ic <- pull.pheno(fake.bc, "sex")

# scan with additive but not the interactive covariate
out.acovar <- scanone(fake.bc, addcovar=ac)

# scan with interactive covariate
out.icovar <- scanone(fake.bc, addcovar=ac, intcovar=ic)

# plot the difference of with and without the interactive covariate
#     This is a LOD score for a test of QTL x covariate interaction
plot(out.icovar-out.acovar)
}

\keyword{arith}
