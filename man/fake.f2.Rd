\name{fake.f2}
\alias{fake.f2}
\docType{data}

\title{Simulated data for an F2 intercross}

\description{ 
  Simulated data for an F2 intercross, obtained using
  \code{\link{sim.cross}}.   
} 

\usage{data(fake.f2)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details. 
}

\details{
There are 200 F2 individuals typed at 94 markers, including 3 on the X
chromosome.  There is one quantitative phenotype, along with an
indication of sex (0=female, 1=male) and the direction of the cross (pgm
= paternal grandmother, 0=A, meaning the cross was (AxB)x(AxB), and
1=B, meaning the cross was (AxB)x(BxA)).

Note that the X chromosome genotypes are coded in a special way (see
\code{\link{read.cross}}).  For the individuals with pgm=0, sex=0,
1=AA and 2=AB; for individuals with pgm=0, sex=1, 1=A and 2=B
(hemizygous); for individuals with pgm=1, sex=0, 1=BB and 2=AB; for
individuals with pgm=1, sex=1, 1=A and 2=B.  \bold{This requires special
care!}

The data were simulated using an additive model with three QTLs on
chromosome 1 (at 30, 50 and 70 cM), one QTL on chromosome 13 (at 30
cM), and one QTL on the X chromosome (at 10 cM). 
}


\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
summary(fake.f2)
plot(fake.f2)
}

\seealso{ \code{\link{sim.cross}}, \code{\link{fake.bc}},
  \code{\link{fake.4way}}, \code{\link{listeria}},
  \code{\link{hyper}},
  \code{\link{bristle3}}, \code{\link{bristleX}} } 

\keyword{datasets}
