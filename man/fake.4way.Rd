\name{fake.4way}
\alias{fake.4way}
\docType{data}

\title{Simulated data for a 4-way cross}

\description{ 
  Simulated data for a phase-known 4-way cross, obtained using
  \code{\link{sim.cross}}.
} 

\usage{data(fake.4way)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details.  
}

\details{
There are 250 individuals typed at 157 markers, including 8 on the X
chromosome.  

There are two phenotypes (including sex, for which 0=female and 1=male).
The quantitative phenotype is affected by three QTLs: two on chromosome
2 at positions 10 and 25 cM on the female genetic map, and one on
chromosome 7 at position 40 cM on the female map. 
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.4way)
\dontshow{fake.4way <- subset(fake.4way,chr=c(2,7))}
plot(fake.4way)
summary(fake.4way)

# estimate recombination fractions
fake.4way <- est.rf(fake.4way)
plotRF(fake.4way)

# estimate genetic maps
ssmap <- est.map(fake.4way, verbose=TRUE)
samap <- est.map(fake.4way, sex.sp=FALSE, verbose=TRUE)
plot(ssmap, samap)

# error lod scores
fake.4way <- calc.genoprob(fake.4way, err=0.01)
fake.4way <- calc.errorlod(fake.4way, err=0.01)
top.errorlod(fake.4way, cutoff=2.5)

# genome scan
fake.4way <- calc.genoprob(fake.4way, step=2.5)
out.hk <- scanone(fake.4way, method="hk")
out.em <- scanone(fake.4way, method="em")
plot(out.em,out.hk,chr=c(2,7))
}

\seealso{ \code{\link{sim.cross}}, \code{\link{fake.bc}},
  \code{\link{fake.f2}}, \code{\link{listeria}},
  \code{\link{hyper}},
  \code{\link{bristle3}}, \code{\link{bristleX}} }

\keyword{datasets}
