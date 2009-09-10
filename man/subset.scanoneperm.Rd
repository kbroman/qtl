\name{subset.scanoneperm}
\alias{subset.scanoneperm}

\title{Subsetting permutation test results}

\description{
  Pull out results for a specified set LOD columns from 
  permutation results from \code{\link[qtl]{scanone}}.
}

\usage{
\method{subset}{scanoneperm}(x, lodcolumn, \dots)
}

\arguments{
 \item{x}{Permutation results from 
   \code{\link[qtl]{scanone}}, run with \code{n.perm}>0.}
 \item{lodcolumn}{A vector specifying which LOD columns to keep (or, if
   negative), omit.  These should be between 1 and the number of LOD
   columns in the input \code{x}.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{scanone} permutation results, but with only the specified
  subset of the data.}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=16:19)}
fake.bc <- calc.genoprob(fake.bc, step=5)
operm <- scanone(fake.bc, method="hk", pheno.col=1:2, n.perm=25)
operm2 <- subset(operm, lodcolumn=2)
}

\seealso{ \code{\link[qtl]{summary.scanoneperm}},
  \code{\link[qtl]{scanone}}, \code{\link[qtl]{c.scanoneperm}},
  \code{\link[qtl]{cbind.scanoneperm}},
  \code{\link[=c.scanoneperm]{rbind.scanoneperm}} }

\keyword{manip}
