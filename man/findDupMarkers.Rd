\name{findDupMarkers}
\alias{findDupMarkers}

\title{Find markers with identical genotype data}

\description{
  Identify sets of markers with identical genotype data.
}

\usage{
findDupMarkers(cross, chr, exact.only=TRUE, adjacent.only=FALSE)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Optional vector specifying which chromosomes to consider.
   This may be a logical, numeric, or character string vector.}
 \item{exact.only}{If TRUE, look only for markers that have matching
 genotypes and the same pattern of missing data; if FALSE, also look for
 cases where one the observed genotypes at marker match those at
 another, and where the first marker has missing genotype whenever the
 genotype for the second marker is missing.}
 \item{adjacent.only}{If TRUE, look only for sets of markers that are
  adjacent to each other.}
}

\value{A list of marker names; each component is a set of markers whose
  genotypes match one other marker, and the name of the component is the
  name of the marker that they match.}


\details{
  If \code{exact.only=TRUE}, we look only for groups of markers whose
  pattern of missing data and observed genotypes match exactly.  One
  marker (chosen at random) is selected as the name of the group (in the
  output of the function).

  If \code{exact.only=FALSE}, we look also for markers whose observed genotypes
  are contained in the observed genotypes of another marker.  We use a
  pair of nested loops, working from the markers with the most observed
  genotypes to the markers with the fewest observed genotypes.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(hyper)

hyper <- drop.nullmarkers(hyper)

dupmar <- findDupMarkers(hyper) # finds 4 pairs
dupmar.adjonly <- findDupMarkers(hyper, adjacent.only=TRUE) # finds 4 pairs

dupmar.nexact <- findDupMarkers(hyper, exact.only=FALSE, adjacent.only=TRUE) # finds 6 pairs

# one might consider dropping the extra markers
totmar(hyper) # 173 markers
hyper <- drop.markers(hyper, unlist(dupmar.adjonly))
totmar(hyper) # 169 markers
}

\seealso{ \code{\link[qtl]{drop.nullmarkers}},
  \code{\link[qtl]{drop.markers}} }

\keyword{manip}
