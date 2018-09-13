\name{subset.scantwoperm}
\alias{subset.scantwoperm}
\alias{[.scantwoperm}

\title{Subsetting two-dimensional permutation test results}

\description{
  Pull out results for a specified set LOD columns from
  permutation results from \code{\link{scantwo}}.
}

\usage{
\method{subset}{scantwoperm}(x, repl, lodcolumn, \dots)
\method{[}{scantwoperm}(x, repl, lodcolumn)
}

\arguments{
 \item{x}{Permutation results from
   \code{\link{scantwo}}, run with \code{n.perm}>0.}
 \item{repl}{A vector specifying which permutation replicates to keep or (if
 negative) omit. Ignored in case of X-chr specific permutations}
 \item{lodcolumn}{A vector specifying which LOD columns to keep or (if
   negative) omit.  These should be between 1 and the number of LOD
   columns in the input \code{x}.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{scantwo} permutation results, but with only the specified
  subset of the data.}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=18:19)}
fake.bc <- calc.genoprob(fake.bc, step=0)
operm <- scantwo(fake.bc, method="hk", pheno.col=1:2, n.perm=5)
operm2 <- subset(operm, lodcolumn=2)

# alternatively
operm2alt <- operm[,2]
}

\seealso{ \code{\link{summary.scantwoperm}},
  \code{\link{scantwo}}, \code{\link{c.scantwoperm}},
  \code{\link[=c.scantwoperm]{rbind.scantwoperm}} }

\keyword{manip}
