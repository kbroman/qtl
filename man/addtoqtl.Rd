\name{addtoqtl}
\alias{addtoqtl}

\title{Add to a qtl object}

\description{
  Add a QTL or multiple QTL to a qtl object.
}

\usage{
addtoqtl(cross, qtl, chr, pos, qtl.name, drop.lod.profile=TRUE)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{qtl}{The qtl object to which additional QTL are to be added.}
 \item{chr}{Vector indicating the chromosome for each new QTL. (These
 should be character strings referring to the chromosomes by name.)}
 \item{pos}{Vector (of same length as \code{chr}) indicating the
   positions on the chromosome for each new QTL. If there is no marker or
   pseudomarker at a position, the nearest position is used.}
 \item{qtl.name}{Optional user-specified name for each new QTL, used in the
   drop-one-term ANOVA table in \code{\link{fitqtl}}.
   If unspecified, the names will be of the form \code{"Chr1@10"} for a
   QTL on Chromsome 1 at 10 cM.} 
 \item{drop.lod.profile}{If TRUE, remove any LOD profiles from the object.}
}

\value{
  An object of class \code{qtl}, just like the input \code{qtl} object,
  but with additional QTL added.  See \code{\link{makeqtl}} for
  details. 
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(qc, 14))}
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")
qtl <- addtoqtl(fake.f2, qtl, 14, 35)
}

\seealso{ \code{\link{makeqtl}}, \code{\link{fitqtl}},
  \code{\link{dropfromqtl}}, \code{\link{replaceqtl}},
  \code{\link{reorderqtl}}  }

\keyword{utilities}
