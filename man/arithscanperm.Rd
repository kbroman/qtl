\name{arithscanperm}
\alias{+.scanoneperm}
\alias{-.scanoneperm}
\alias{+.scantwoperm}
\alias{-.scantwoperm}

\title{Arithmetic Operators for permutation results}

\description{
  Add or subtract LOD scores in permutation results from
  \code{\link{scanone}} or \code{\link{scantwo}}.
}

\usage{
\special{perm1+perm2}
\special{perm1-perm2}
}

\arguments{
  \item{perm1, perm2}{Permutation results from
    \code{\link{scanone}} or \code{\link{scantwo}}, on the
    same set of chromosomes and markers.}
}

\details{
  This is used to calculate the sum or difference of LOD scores of two
  sets of permutation results from \code{\link{scanone}} or
  \code{\link{scantwo}}.  One must be careful to ensure that the
  permutations are perfectly linked, which
  will require the use of \code{\link[base:Random]{set.seed}}.
}

\value{
  The same data structure as the input objects, with LOD scores
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

# set seed
theseed <- round(runif(1, 1, 10^8))
set.seed(theseed)

# permutations with additive but not the interactive covariate
\dontrun{operm.acovar <- scanone(fake.bc, addcovar=ac, n.perm=1000)
}\dontshow{operm.acovar <- scanone(fake.bc, pheno.col=1, method="hk", addcovar=ac, n.perm=100)
}

# re-set the seed
set.seed(theseed)

# permutations with interactive covariate
\dontrun{operm.icovar <- scanone(fake.bc, addcovar=ac, intcovar=ic,
                      n.perm=1000)
}\dontshow{operm.icovar <- scanone(fake.bc, pheno.col=1, addcovar=ac,
                      intcovar=ic, method="hk", n.perm=100)
}

# permutation results for the QTL x covariate interaction
operm.gxc <- operm.icovar - operm.acovar

# LOD thresholds
summary(operm.gxc)
}

\keyword{arith}
