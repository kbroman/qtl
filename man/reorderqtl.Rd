\name{reorderqtl}
\alias{reorderqtl}

\title{Reorder the QTL in a qtl object}

\description{
  This function changes the order of the QTL in a QTL object.
}

\usage{
reorderqtl(qtl, neworder)
}

\arguments{
 \item{qtl}{A qtl object, as created by \code{\link{makeqtl}}.}
 \item{neworder}{A vector containing the positive integers up to the
   number of QTL in \code{qtl}, indicating the new order for the QTL.
   If missing, the QTL are ordered by chromosome and then by their
   position within a chromosome.}
}


\value{
  The input \code{qtl} object, with the loci reordered.
}

\details{
  Everything in the input \code{qtl} is reordered except the
  \code{altname} component, which contains names of the form \code{Q1},
  \code{Q2}, etc.  
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

qtl <- reorderqtl(qtl, c(2,3,1))
qtl

qtl <- reorderqtl(qtl)
qtl
}

\seealso{ \code{\link{makeqtl}}, \code{\link{fitqtl}},
  \code{\link{dropfromqtl}}, \code{\link{addtoqtl}},
  \code{\link{replaceqtl}} }

\keyword{utilities}
