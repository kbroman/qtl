\name{locateXO}
\alias{locateXO}

\title{Estimate locations of crossovers}

\description{
  Estimate the locations of crossovers for each individual on a given
  chromosome.
}

\usage{
locateXO(cross, chr, full.info=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Chromosome to investigate (if unspecified, the first
   chromosome is considered).
   This should be a character string referring to a chromosome
   by name; numeric values are converted to strings.}
 \item{full.info}{If TRUE, output will include information on the left
   and right endpoints of the intervals to which recombination events
   are known, as well as the corresponding marker indices.}
}

\details{
  For each individual we detemine the locations of obligate crossovers,
  and estimate their location to be at the midpoint between the nearest
  flanking typed markers.
  
  The function currently only works for a backcross, intercross, or
  recombinant inbred line.
}

\value{
  A list with one component per individual.  Each component is either
  NULL or is a numeric vector with the estimated crossover locations.

  If \code{full.info=TRUE}, in place of a numeric vector with estimated
  locations, there is a matrix that includes those locations, the left
  and right endpoints of the intervals to which crossovers can be
  placed, the marker indices corresponding to those endpoint, and
  genotype codes for the genotypes to the left and right of each crossover.  The
  final column indicates the number of typed markers between the current
  crossover and the next one (useful for identifying potential
  genotyping errors).
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(hyper)
xoloc <- locateXO(hyper, chr=4)
table(sapply(xoloc, length))
}

\seealso{ \code{\link{countXO}}, \code{\link{cleanGeno}} }

\keyword{utilities}
