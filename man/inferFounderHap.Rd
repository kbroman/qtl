\name{inferFounderHap}
\alias{inferFounderHap}

\title{Crude reconstruction of founder haplotypes in multi-parent RIL}

\description{
  Uses groups of adjacent markers to infer the founder haplotypes in SNP
  data on multi-parent recombinant inbred lines.
}

\usage{
inferFounderHap(cross, chr, max.n.marker=15)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Indicator of chromosome to consider.  If multiple
 chromosomes are selected, only the first is used.}
 \item{max.n.marker}{Maximum number of adjacent markers to consider.}
}

\details{
}

\value{
  A matrix of dimension \code{nind(cross)} \eqn{\times}{x} no. markers,
  with the inferred founder origin for each line at each marker. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
map <- sim.map(100, n.mar=101, include.x=FALSE, eq.spacing=TRUE)
founderGeno <- simFounderSnps(map, "8")
ril <- sim.cross(map, n.ind=10, type="ri8sib", founderGeno=founderGeno)

h <- inferFounderHap(ril, max.n.marker=11)
mean(!is.na(h)) # proportion inferred
plot(map[[1]], h[1,], ylim=c(0.5, 8.5), xlab="Position", ylab="Genotype")
}

\seealso{ \code{\link{sim.geno}}, \code{\link{calc.genoprob}},
          \code{\link{fill.geno}}, \code{\link{argmax.geno}} }

\keyword{utilities}
