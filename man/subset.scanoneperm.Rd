\name{subset.scanoneperm}
\alias{subset.scanoneperm}
\alias{[.scanoneperm}

\title{Subsetting permutation test results}

\description{
  Pull out results for a specified set LOD columns from
  permutation results from \code{\link{scanone}}.
}

\usage{
\method{subset}{scanoneperm}(x, repl, lodcolumn, \dots)
\method{[}{scanoneperm}(x, repl, lodcolumn)
}

\arguments{
 \item{x}{Permutation results from
   \code{\link{scanone}}, run with \code{n.perm}>0.}
 \item{repl}{A vector specifying which permutation replicates to keep or (if
 negative) omit.}
 \item{lodcolumn}{A vector specifying which LOD columns to keep or (if
   negative) omit.  These should be between 1 and the number of LOD
   columns in the input \code{x}.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{scanone} permutation results, but with only the specified
  subset of the data.}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=16:19)}
fake.bc <- calc.genoprob(fake.bc, step=5)
operm <- scanone(fake.bc, method="hk", pheno.col=1:2, n.perm=25)
operm2 <- subset(operm, lodcolumn=2)

# alternatively
operm2alt <- operm[,2]
}

\seealso{ \code{\link{summary.scanoneperm}},
  \code{\link{scanone}}, \code{\link{c.scanoneperm}},
  \code{\link{cbind.scanoneperm}},
  \code{\link[=c.scanoneperm]{rbind.scanoneperm}} }

\keyword{manip}
