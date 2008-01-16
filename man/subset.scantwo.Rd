\name{subset.scantwo}
\alias{subset.scantwo}

\title{Subsetting the results of a 2-d genome scan}

\description{
  Pull out a specified set of chromosomes and/or LOD columns from 
  \code{\link[qtl]{scantwo}} output.
}

\usage{
\method{subset}{scantwo}(x, chr, lodcolumn, \dots)
}

\arguments{
 \item{x}{An object of class \code{scantwo}, output from
   \code{\link[qtl]{scantwo}}.}
 \item{chr}{A character vector specifying which chromosomes to keep.}
 \item{lodcolumn}{A vector specifying which LOD columns to keep (or, if
   negative), omit.  These should be between 1 and the number of LOD
   columns in the input \code{x}.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{scantwo} object, but with only the specified
  subset of the data.}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=16:19)}
fake.bc <- calc.genoprob(fake.bc)
out <- scantwo(fake.bc, method="hk", pheno.col=1:2)

summary(subset(out, chr=18:19))
}

\seealso{ \code{\link[qtl]{summary.scantwo}},
  \code{\link[qtl]{scantwo}} }

\keyword{manip}
