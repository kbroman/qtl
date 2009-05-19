\name{clean.cross}
\alias{clean.cross}

\title{Remove derived data}

\description{
  Remove any intermediate calculations from a cross object.
}

\usage{
\method{clean}{cross}(object, \dots)
}

\arguments{
 \item{object}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{\dots}{Ignored at this point.}
}

\value{
  The input object, with any intermediate calculations
  (such as is produced by \code{\link[qtl]{calc.genoprob}},
  \code{\link[qtl]{argmax.geno}} and \code{\link[qtl]{sim.geno}})
  removed.
}

\examples{
data(fake.f2)
names(fake.f2$geno)
fake.f2 <- calc.genoprob(fake.f2)
names(fake.f2$geno)
fake.f2 <- clean(fake.f2)
names(fake.f2$geno)
}

\seealso{ \code{\link[qtl]{drop.nullmarkers}},
  \code{\link[qtl]{drop.markers}}, \code{\link[qtl]{clean.scantwo}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{manip}
