\name{subset.scantwo}
\alias{subset.scantwo}

\title{Subsetting the results of a 2-d genome scan}

\description{
  Pull out a specified set of chromosomes and/or LOD columns from 
  \code{\link{scantwo}} output.
}

\usage{
\method{subset}{scantwo}(x, chr, lodcolumn, \dots)
}

\arguments{
 \item{x}{An object of class \code{scantwo}, output from
   \code{\link{scantwo}}.}
 \item{chr}{Optional vector specifying which chromosomes to keep.
 This should be a vector of character strings referring to
 chromosomes by name; numeric values are converted to strings.  Refer to
 chromosomes with a preceding \code{-} to have all chromosomes but those
 considered.  A logical (TRUE/FALSE) vector may also be used.} 
 \item{lodcolumn}{A vector specifying which LOD columns to keep (or, if
   negative), omit.  These should be between 1 and the number of LOD
   columns in the input \code{x}.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{scantwo} object, but with only the specified
  subset of the data.}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=16:19)}
fake.bc <- calc.genoprob(fake.bc)
out <- scantwo(fake.bc, method="hk", pheno.col=1:2)

summary(subset(out, chr=18:19))
}

\seealso{ \code{\link{summary.scantwo}},
  \code{\link{scantwo}} }

\keyword{manip}
