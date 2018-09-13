\name{summary.fitqtl}
\alias{summary.fitqtl}

\title{Summary of fit of qtl model}

\description{
  Print summary information about the results of \code{\link{fitqtl}}.
}

\usage{
\method{summary}{fitqtl}(object, pvalues=TRUE, simple=FALSE, \dots)
}
\arguments{
 \item{object}{Output from \code{\link{fitqtl}}.}
 \item{pvalues}{If FALSE, don't include p-values in the summary.}
 \item{simple}{If TRUE, don't include p-values or sums of squares in the summary.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.fitqtl}, which is not all that
  different than the input, but when printed gives summary information
  about the results.
}

\author{Hao Wu; Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 8, 13)
qp <- c(26, 56, 28)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

# fit model with 3 interacting QTLs interacting
# (performing a drop-one-term analysis)
lod <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1*Q2*Q3,
              method="hk")
summary(lod)
}

\seealso{ \code{\link{fitqtl}}, \code{\link{makeqtl}},
\code{\link{scanqtl}} }  

\keyword{print}
