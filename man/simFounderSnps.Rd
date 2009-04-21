\name{simFounderSnps}
\alias{simFounderSnps}

\title{Simulate founder SNPs for a multiple-strain RIL}

\description{
  Simulate genotype data for the founding strains for a panel of
  multiple-strain RIL.
}

\usage{
simFounderSnps(map, n.str=c("4","8"), pat.freq)
}

\arguments{
 \item{map}{A list whose components are vectors containing the marker
   locations on each of the chromosomes.}
 \item{n.str}{Number of founding strains (4 or 8).}
 \item{pat.freq}{Frequency of SNP genotype patterns in the founder (a
   vector of length \code{n.str}/2 + 1): (monoallelic, SNP unique to one
   founder, SNP present in 2 founders, [and, for the case of 8 founders,
   SNP in 3/8 founders, SNP in 4/8 founders].)}
}

\details{
  The SNPs are simulated to be in linkage equilibrium. 
}  

\value{
  A vector of the same length as there are chromosomes in \code{map},
  with each component being a matrix of 0's and 1's, of dim \code{n.str}
  x \code{n.mar}. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\seealso{ \code{\link[qtl]{sim.map}}, \code{\link[qtl]{sim.cross}} }

\examples{
data(map10)
x <- simFounderSnps(map10, "8", c(0, 0.5, 0.2, 0.2, 0.1))
}
                    
\keyword{datagen}
