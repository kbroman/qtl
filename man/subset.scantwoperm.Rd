\name{subset.scantwoperm}
\alias{subset.scantwoperm}

\title{Subsetting two-dimensional permutation test results}

\description{
  Pull out results for a specified set LOD columns from 
  permutation results from \code{\link[qtl]{scantwo}}.
}

\usage{
\method{subset}{scantwoperm}(x, lodcolumn, \dots)
}

\arguments{
 \item{x}{Permutation results from 
   \code{\link[qtl]{scantwo}}, run with \code{n.perm}>0.}
 \item{lodcolumn}{A vector specifying which LOD columns to keep (or, if
   negative), omit.  These should be between 1 and the number of LOD
   columns in the input \code{x}.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{scantwo} permutation results, but with only the specified
  subset of the data.}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=18:19)}
fake.bc <- calc.genoprob(fake.bc, step=0)
operm <- scantwo(fake.bc, method="hk", pheno.col=1:2, n.perm=5)
operm2 <- subset(operm, lodcolumn=2)
}

\seealso{ \code{\link[qtl]{summary.scantwoperm}},
  \code{\link[qtl]{scantwo}}, \code{\link[qtl]{c.scantwoperm}},
  \code{\link[qtl]{cbind.scantwoperm}},
  \code{\link[qtl]{rbind.scantwoperm}} }

\keyword{manip}
