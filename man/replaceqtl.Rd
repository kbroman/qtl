\name{replaceqtl}
\alias{replaceqtl}

\title{Replace a QTL in a qtl object with a different position}

\description{
  This function replaces a QTL or QTLs in a qtl object with a different
  position.
}

\usage{
replaceqtl(cross, qtl, index, chr, pos, qtl.name, drop.lod.profile=TRUE)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{qtl}{A qtl object, as created by \code{\link{makeqtl}}.}
 \item{index}{Numeric index indicating the QTL to be replaced.}
 \item{chr}{Vector (of same length as \code{index}) indicating the
   chromosomes for the new QTL.}
 \item{pos}{Vector (of same length as \code{index}) indicating the
   positions for the new QTL. If there is no marker or pseudomarker at a
 position, the nearest position is used.}
 \item{qtl.name}{Optional vector (of same length as \code{index})
   of user-specified names for each new QTL, used in the
   drop-one-term ANOVA table in \code{\link{fitqtl}}.
   If unspecified, the names will be of the form \code{"Chr1@10"} for a
   QTL on Chromsome 1 at 10 cM.}
 \item{drop.lod.profile}{If TRUE, remove any LOD profiles from the object.}
}


\value{
  The input \code{qtl} object,
  but with some QTL replaced by new ones.  See \code{\link{makeqtl}} for
  details on the format.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

qtl <- replaceqtl(fake.f2, qtl, 2, 6, 48.1)
}

\seealso{ \code{\link{makeqtl}}, \code{\link{fitqtl}},
  \code{\link{dropfromqtl}}, \code{\link{addtoqtl}},
  \code{\link{reorderqtl}} }

\keyword{utilities}
