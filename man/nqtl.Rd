\name{nqtl}
\alias{nqtl}

\title{Determine the number of QTL in a QTL object}

\description{
  Determine the number of QTL in a QTL object.
}

\usage{
nqtl(qtl)
}

\arguments{
 \item{qtl}{An object of class \code{qtl}.  See
   \code{\link{makeqtl}} for details.}
}

\value{
  The number of QTL in the input QTL object.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c("1", "6", "13")
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

nqtl(qtl)
}

\seealso{ \code{\link{makeqtl}}, \code{\link{fitqtl}},
  \code{\link{dropfromqtl}}, \code{\link{replaceqtl}},
  \code{\link{addtoqtl}}, \code{\link{summary.qtl}},
  \code{\link{reorderqtl}} }

\keyword{print}
