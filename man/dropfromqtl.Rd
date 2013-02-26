\name{dropfromqtl}
\alias{dropfromqtl}

\title{Drop a QTL from a qtl object}

\description{
  Drop a QTL or multiple QTL from a QTL object
}

\usage{
dropfromqtl(qtl, index, chr, pos, qtl.name, drop.lod.profile=TRUE)
}

\arguments{
 \item{qtl}{A qtl object, as created by \code{\link{makeqtl}}.}
 \item{index}{Vector specifying the numeric indices of the QTL to be dropped.}
 \item{chr}{Vector indicating the chromosome for each QTL to drop.}
 \item{pos}{Vector (of same length as \code{chr}) indicating the
   positions of the QTL to be dropped.}
 \item{qtl.name}{Vector specifying the names of the QTL to be dropped.}
 \item{drop.lod.profile}{If TRUE, remove any LOD profiles from the object.}
}

\details{
  Provide either \code{chr} and \code{pos}, or one of \code{qtl.name} or
  \code{index}. 
}

\value{
  The input \code{qtl} object with the specified QTL omitted.   See \code{\link{makeqtl}} for
  details on the format. 
}

\author{Karl W Broman, \code{kbroman@biostat.wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

newqtl <- dropfromqtl(qtl, chr=1, pos=25.8)
altqtl <- dropfromqtl(qtl, index=1)
}

\seealso{ \code{\link{makeqtl}}, \code{\link{fitqtl}},
  \code{\link{addtoqtl}}, \code{\link{replaceqtl}} ,
  \code{\link{reorderqtl}}}

\keyword{utilities}
