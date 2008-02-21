\name{countXO}
\alias{countXO}

\title{Count number of obligate crossovers for each individual}

\description{
  Count the number of obligate crossovers for each individual in a
  cross, either by chromosome or overall.
}

\usage{
countXO(cross, chr, bychr=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{The chromosomes to investigate.  If missing, all chromosomes
   are considered.}
 \item{bychr}{If TRUE, return counts for each individual chromosome; if
   FALSE, return the overall number across the selected chromosomes.}
}

\details{
  For each individual we count the minimal number of crossovers that
  explain the observed genotype data.  
}

\value{
  If \code{bychr=TRUE}, a matrix of counts is returned, with rows
  corresponding to individuals and columns corresponding to
  chromosomes. 

  If \code{bychr=FALSE}, a vector of counts (the total number of
  crossovers across all selected chromosomes) is returned.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(hyper)
plot(countXO(hyper))
}

\seealso{ \code{\link[qtl]{ripple}} }

\keyword{utilities}
