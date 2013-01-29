\name{c.scantwoperm}
\alias{c.scantwoperm}
\alias{rbind.scantwoperm}

\title{Combine data from scantwo permutations}

\description{
  Concatenate the data for multiple runs of \code{\link{scantwo}}
  with \code{n.perm > 0}.
}

\usage{
\method{c}{scantwoperm}(\dots)
\method{rbind}{scantwoperm}(\dots)
}

\arguments{
 \item{\dots}{A set of objects of class \code{scantwoperm}. (This can
 also be a list of \code{scantwoperm} objects.)  These are
   the permutation results from \code{\link{scantwo}} (that is,
   when \code{n.perm > 0}).  These must all concern the same number of
   LOD columns. (That is, they must have been created with the same number
   of phenotypes, and it is assumed that they were generated in
   precisely the same way.)}
}

\value{The concatenated input, as a \code{scantwoperm} object.}

\details{
  The aim of this function is to concatenate the results from multiple
  runs of a permutation test \code{\link{scantwo}}, to assist with
  the case that such permutations are done on multiple processors in
  parallel. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"), ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2)
\dontrun{operm1 <- scantwo(fake.f2, method="hk", n.perm=50)
operm2 <- scantwo(fake.f2, method="hk", n.perm=50)}\dontshow{
operm1 <- scantwo(fake.f2, method="hk", n.perm=3)
operm2 <- scantwo(fake.f2, method="hk", n.perm=3)}

operm <- c(operm1, operm2)
}

\seealso{ \code{\link{summary.scantwoperm}},
  \code{\link{scantwo}}, \code{\link{summary.scantwo}},
  \code{\link{c.scanoneperm}} }

\keyword{manip}
