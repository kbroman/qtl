\name{replacemap.scanone}
\alias{replacemap.scanone}

\title{Replace the genetic map in QTL mapping results with an alternate map}

\description{
  Replace the positions of LOD scores in output from
  \code{\link[qtl]{scanone}} with values
  based on an alternative map (such as a physical map), with
  pseudomarker locations determined by linear interpolation.
}

\usage{
\method{replacemap}{scanone}(object, map)
}

\arguments{
 \item{object}{An object of class \code{"scanone"}, as output by the
 function \code{\link[qtl]{scanone}}.}
 \item{map}{A list containing the alternative genetic map.  All
 chromosomes in \code{object} should have corresponding chromosomes in
 \code{map}, and markers must be in the same order in the two maps.
 There must be at least two markers on each chromosome in \code{map}
 that appear in \code{object}.}
}

\value{The input \code{object} with the positions of LOD scores 
  revised to match those in the input \code{map}.}

\details{
  The positions of pseudomarkers are determined by linear interpolation
  between markers.  In the case of pseudomarkers beyond the ends of the
  terminal markers on chromosomes, we use the overall lengths of the
  chromosome in \code{object} and \code{map} to determine the new
  spacing.
}  

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
origmap <- pull.map(fake.f2)
newmap <- est.map(fake.f2)
fake.f2 <- replacemap(fake.f2, newmap)
out <- scanone(fake.f2, method="hk")
out.rev <- replacemap(out, origmap)
}

\seealso{ \code{\link[qtl]{replacemap.cross}},
  \code{\link[qtl]{est.map}}, \code{\link[qtl]{replacemap.scantwo}} }

\keyword{manip}
