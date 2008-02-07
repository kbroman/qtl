\name{summary.fitqtl}
\alias{summary.fitqtl}

\title{Summary of fit of qtl model}

\description{
  Print summary information about the results of \code{\link[qtl]{fitqtl}}.
}

\usage{
\method{summary}{fitqtl}(object, \dots)
}
\arguments{
 \item{object}{Output from \code{\link[qtl]{fitqtl}}.}
 \item{\dots}{Ignored at this point.}
}

%\details{
%  To be added.
%}

\value{
  An object of class \code{summary.fitqtl}, which is not all that
  different than the input, but when printed gives summary information
  about the results.
}

\references{
  Sen, \'S. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Hao Wu}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 8, 13)
qp <- c(26, 56, 28)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp)

# fit model with 3 interacting QTLs interacting
# (performing a drop-one-term analysis)
lod <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1*Q2*Q3)
summary(lod)
}

\seealso{ \code{\link[qtl]{fitqtl}}, \code{\link[qtl]{makeqtl}},
\code{\link[qtl]{scanqtl}} }  

\keyword{print}
