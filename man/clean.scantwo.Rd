\name{clean.scantwo}
\alias{clean.scantwo}

\title{Clean up scantwo output}

\description{
  In an object output from \code{\link[qtl]{scantwo}}, replaces negative
  and missing LOD scores with 0, and replaces LOD scores for pairs of
  positions that are not separated by \code{n.mar} markers, or that are
  less than \code{distance} cM apart, with 0.  Further, if the LOD
  for full model is less than the LOD for the additive model, the
  additive LOD is pasted over the full LOD.
}

\usage{
\method{clean}{scantwo}(object, n.mar=1, distance=0, \dots)
}

\arguments{
 \item{object}{An object of class \code{scantwo}.  See
   \code{\link[qtl]{scantwo}} for details.}
 \item{n.mar}{}
 \item{distance}{}
 \item{\dots}{Ignored at this point.}
}

\value{
  The input scantwo object, with any negative or missing LOD scores
  replaced by 0, and LOD scores for pairs of positions separated by
  fewer than \code{n.mar} markers, or less than \code{distance} cM, are
  set to 0.
  Also, if the LOD for the full model is less than the LOD for the
  additive model, the additive LOD is used in place of the full LOD.
}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
fake.f2 <- calc.genoprob(fake.f2, step=5)
out2 <- scantwo(fake.f2, method="hk")
out2 <- clean(out2)
out2cl2 <- clean(out2, n.mar=2, distance=5)
}

\seealso{ \code{\link[qtl]{scantwo}},
  \code{\link[qtl]{summary.scantwo}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{manip}
