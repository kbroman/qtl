\name{inferFounderHap}
\alias{inferFounderHap}

\title{Crude reconstruction of founder haplotypes in multi-parent RIL}

\description{
  Uses groups of adjacent markers to infer the founder haplotypes in SNP
  data on multi-parent recombinant inbred lines.
}

\usage{
inferFounderHap(cross, chr, max.n.markers=15)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Indicator of chromosome to consider.  If multiple
 chromosomes are selected, only the first is used.}
 \item{max.n.markers}{Maximum number of adjacent markers to consider.}
}

\details{
  We omit SNPs for which any of the founders are missing.

  We then consider groups of adjacent SNPs, looking for founder
  haplotypes that are unique; RIL sharing such a unique haplotype are
  then inferred to have that founder's DNA.

  We consider each marker as the center of a haplotype, and consider
  haplotypes of size 1, 3, 5, \dots, \code{max.n.markers}.  We end the
  extension of the haplotypes when all founders have a unique haplotype.
}

\value{
  A matrix of dimension \code{nind(cross)} \eqn{\times}{x} no. markers,
  with the inferred founder origin for each line at each marker.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
map <- sim.map(100, n.mar=101, include.x=FALSE, eq.spacing=TRUE)
founderGeno <- simFounderSnps(map, "8")
ril <- sim.cross(map, n.ind=10, type="ri8sib", founderGeno=founderGeno)

h <- inferFounderHap(ril, max.n.markers=11)
mean(!is.na(h)) # proportion inferred
plot(map[[1]], h[1,], ylim=c(0.5, 8.5), xlab="Position", ylab="Genotype")
}

\seealso{ \code{\link{sim.geno}}, \code{\link{calc.genoprob}},
          \code{\link{fill.geno}}, \code{\link{argmax.geno}} }

\keyword{utilities}
