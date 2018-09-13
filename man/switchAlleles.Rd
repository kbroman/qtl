\name{switchAlleles}
\alias{switchAlleles}

\title{Switch alleles at selected markers}

\description{
  Switch alleles at selected markers in a cross object.
}

\usage{switchAlleles(cross, markers, switch=c("AB", "CD", "ABCD", "parents"))}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{markers}{Names of markers whose alleles are to be switched.}
 \item{switch}{For a 4-way cross, indicates how to switch the alleles (A
 for B, C for D, both A for B and C for D), or both A for C and B for D (\code{parents}).}
}

\value{
  The input cross object, with alleles at selected markers switched.
}

\details{
  For a backcross, we exchange homozygotes (AA) and heterozygotes (AB).

  For doubled haploids and recombinant inbred lines, we
  exchange the two homozygotes.

  For an intercross, we exchange the two homozygotes, and exchange C
  (i.e., not AA) and D (i.e., not BB). (The heterozygotes in an intercross are
  left unchanged.)

  For a 4-way cross, we consider the argument \code{switch}, and the
  exchanges among the genotypes are more complicated.
}


\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
geno.crosstab(fake.f2, "D5M391", "D5M81")

# switch homozygotes at marker D5M391
fake.f2 <- switchAlleles(fake.f2, "D5M391")

geno.crosstab(fake.f2, "D5M391", "D5M81")

\dontrun{fake.f2 <- est.rf(fake.f2)
checkAlleles(fake.f2)
}
}

\seealso{ \code{\link{checkAlleles}}, \code{\link{est.rf}}, \code{\link{geno.crosstab}} }

\keyword{utilities}
