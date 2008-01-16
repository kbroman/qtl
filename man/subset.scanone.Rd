\name{subset.scanone}
\alias{subset.scanone}

\title{Subsetting the results of a genome scan}

\description{
  Pull out a specified set of chromosomes and/or LOD columns from 
  \code{\link[qtl]{scanone}} output.
}

\usage{
\method{subset}{scanone}(x, chr, lodcolumn, \dots)
}

\arguments{
 \item{x}{An object of class \code{scanone}, output from
   \code{\link[qtl]{scanone}}.}
 \item{chr}{A character vector specifying which chromosomes to keep.}
 \item{lodcolumn}{A vector specifying which LOD columns to keep (or, if
   negative), omit.  These should be between 1 and the number of LOD
   columns in the input \code{x}.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{scanone} object, but with only the specified
  subset of the data.}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=16:19)}
fake.bc <- calc.genoprob(fake.bc, step=2.5)
out <- scanone(fake.bc, method="hk", pheno.col=1:2)

summary(subset(out, chr=18:19), format="allpeaks")
}

\seealso{ \code{\link[qtl]{summary.scanone}},
  \code{\link[qtl]{scanone}} }

\keyword{manip}
