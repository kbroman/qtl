\name{addtoqtl}
\alias{addtoqtl}

\title{Add to a qtl object}

\description{
  Add a QTL or multiple QTL to a qtl object.
}

\usage{
addtoqtl(cross, qtl, chr, pos, qtl.name)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{qtl}{The qtl object to which additional QTL are to be added.}
 \item{chr}{Vector indicating the chromosome for each new QTL.}
 \item{pos}{Vector (of same length as \code{chr}) indicating the
   positions on the chromosome for each new QTL. If there is no marker or
   pseudomarker at a position, the nearest position is used.}
 \item{qtl.name}{Optional user-specified name for each new QTL, used in the
   drop-one-term ANOVA table in \code{\link[qtl]{fitqtl}}.
   If unspecified, the names will be of the form \code{"Chr1@10"} for a
   QTL on Chromsome 1 at 10 cM.} 
}

\value{
  An object of class \code{qtl}, just like the input \code{qtl} object,
  but with additional QTL added.  See \code{\link[qtl]{makeqtl}} for
  details. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(qc, 14))}
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="draws")
qtl <- addtoqtl(fake.f2, qtl, 14, 35)
}

\seealso{ \code{\link[qtl]{makeqtl}}, \code{\link[qtl]{fitqtl}},
  \code{\link[qtl]{dropfromqtl}}, \code{\link[qtl]{replaceqtl}},
  \code{\link[qtl]{reorderqtl}}  }

\keyword{utilities}
