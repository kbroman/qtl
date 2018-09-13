\name{fake.bc}
\alias{fake.bc}
\docType{data}

\title{Simulated data for a backcross}

\description{
  Simulated data for a backcross, obtained using
  \code{\link{sim.cross}}.
}

\usage{data(fake.bc)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details.
}


\details{
  There are 400 backcross individuals typed at 91 markers and with two
  phenotypes and two covariates (sex and age).

  The two phenotypes are due to four QTLs, with no epistasis.  There is
  one on chromosome 2 (at 30 cM), two on chromosome 5 (at 10 and 50 cM),
  and one on chromosome 10 (at 30 cM).  The QTL on chromosome 2 has an
  effect only in the males (sex=1); the two QTLs on chromosome 5 have
  effect in coupling for the first phenotype and in repulsion for the
  second phenotype. Age has an effect of increasing the phenotypes.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc,chr=c(2,5,10))}
summary(fake.bc)
plot(fake.bc)

# genome scans without covariates
fake.bc <- calc.genoprob(fake.bc, step=2.5)
out.nocovar <- scanone(fake.bc, pheno.col=1:2)

# genome scans with covariates
ac <- pull.pheno(fake.bc, c("sex","age"))
ic <- pull.pheno(fake.bc, "sex")
out.covar <- scanone(fake.bc, pheno.col=1:2,
                      addcovar=ac, intcovar=ic)

# summaries
summary(out.nocovar, thr=3, format="allpeaks")
summary(out.covar, thr=3, format="allpeaks")

# plots
plot(out.nocovar, out.covar, chr=c(2,5,10), lod=1, col="blue",
     lty=1:2, ylim=c(0,13))
plot(out.nocovar, out.covar, chr=c(2,5,10), lod=2, col="red",
     lty=1:2, add=TRUE)
}

\seealso{ \code{\link{sim.cross}}, \code{\link{fake.4way}},
  \code{\link{fake.f2}}, \code{\link{listeria}},
  \code{\link{hyper}},
  \code{\link{bristle3}}, \code{\link{bristleX}} }

\keyword{datasets}
