\name{summary.qtl}
\alias{summary.qtl}

\title{Print summary of a QTL object}

\description{
  Print summary information about a \code{qtl} object.
}

\usage{
\method{summary}{qtl}(object, \dots)
}

\arguments{
 \item{object}{An object of class \code{qtl}, created by \code{\link[qtl]{makeqtl}}.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.qtl}, which is just a data.frame 
  containing the chromosomes, positions, and number of possible
  genotypes for each QTL.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 6, 13)
qp <- c(25.8, 33.6, 18.63)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="draws")

summary(qtl)
}

\seealso{ \code{\link[qtl]{makeqtl}} } 

\keyword{print}
