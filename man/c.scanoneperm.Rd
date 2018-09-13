\name{c.scanoneperm}
\alias{c.scanoneperm}
\alias{rbind.scanoneperm}

\title{Combine data from scanone permutations}

\description{
  Concatenate the data for multiple runs of \code{\link{scanone}}
  with \code{n.perm > 0}.
}

\usage{
\method{c}{scanoneperm}(\dots)
\method{rbind}{scanoneperm}(\dots)
}

\arguments{
 \item{\dots}{A set of objects of class \code{scanoneperm}. (This can
 also be a list of \code{scanoneperm} objects.) These are
   the permutation results from \code{\link{scanone}} (that is,
   when \code{n.perm > 0}).  These must all have the same number of
   columns. (That is, they must have been created with the same number
   of phenotypes, and it is assumed that they were generated in
   precisely the same way.)}
}

\value{The concatenated input, as a \code{scanoneperm} object.}

\details{
  The aim of this function is to concatenate the results from multiple
  runs of a permutation test \code{\link{scanone}}, to assist with
  the case that such permutations are done on multiple processors in
  parallel. 
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"))}
fake.f2 <- calc.genoprob(fake.f2)
operm1 <- scanone(fake.f2, method="hk", n.perm=100, perm.Xsp=TRUE)
operm2 <- scanone(fake.f2, method="hk", n.perm=50, perm.Xsp=TRUE)

operm <- c(operm1, operm2)
}

\seealso{ \code{\link{summary.scanoneperm}},
  \code{\link{scanone}}, \code{\link{cbind.scanoneperm}},
  \code{\link{c.scantwoperm}} }

\keyword{manip}
