\name{cbind.scantwoperm}
\alias{cbind.scantwoperm}

\title{Combine scantwo permutations by column}

\description{
  Column-bind permutations results from \code{\link{scantwo}}
  for multiple phenotypes or models.
}

\usage{
\method{cbind}{scantwoperm}(\dots)
}

\arguments{
 \item{\dots}{A set of objects of class \code{scantwoperm}. (This can
   also be a list of \code{scantwoperm} objects.)  These are
   the permutation results from \code{\link{scantwo}} (that is,
   when \code{n.perm > 0}).  These must all concern the same number of
   permutations.}
}

\value{The column-binded input, as a \code{scantwoperm} object.}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=c(18:19), ind=1:50)}
fake.bc <- calc.genoprob(fake.bc)
\dontrun{operm1 <- scantwo(fake.bc, pheno.col=1, method="hk", n.perm=50)
operm2 <- scantwo(fake.bc, pheno.col=2, method="hk", n.perm=50)}\dontshow{
operm1 <- scantwo(fake.bc, pheno.col=1, method="hk", n.perm=3)
operm2 <- scantwo(fake.bc, pheno.col=2, method="hk", n.perm=3)}

operm <- cbind(operm1, operm2)
}

\seealso{ \code{\link{scantwo}}, \code{\link{c.scantwoperm}},
  \code{\link{summary.scantwoperm}} }

\keyword{manip}
