\name{sim.map}
\alias{sim.map}

\title{Simulate a genetic map}

\description{
  Simulate the positions of markers on a genetic map.
}

\usage{
sim.map(len=rep(100,20), n.mar=10, anchor.tel=TRUE,
        include.x=TRUE, sex.sp=FALSE, eq.spacing=FALSE)
}

\arguments{
 \item{len}{A vector specifying the chromosome lengths (in cM)}
 \item{n.mar}{A vector specifying the number of markers per chromosome.}
 \item{anchor.tel}{If true, markers at the two telomeres will always be
   included, so if \code{n.mar} = 1 or 2, we'll give just the two
   telomeric markers.}
 \item{include.x}{Indicates whether the last chromosome should be
   considered the X chromosome.}
 \item{sex.sp}{Indicates whether to create sex-specific maps, in which
   case the output will be a vector of 2-row matrices, with rows
   corresponding to the maps for the two sexes.}
 \item{eq.spacing}{If TRUE, markers will be equally spaced.}
}

\details{
  Aside from the telomeric markers, marker positions are simulated as
  iid Uniform(\eqn{0,L}).  If \code{len} or \code{n.mar} has just one element,
  it is expanded to the length of the other argument.  If they both have
  just one element, only one chromosome is simulated.

  If \code{eq.spacing} is TRUE, markers are equally spaced between 0 and
  \eqn{L}.  If \code{anchor.tel} is FALSE, telomeric markers are not
  included.
}

\value{
  A list of vectors, each specifying the locations of the markers.  Each
  component of the list is given class \code{A} or \code{X}, according
  to whether it is autosomal or the X chromosome.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
# simulate 4 autosomes, each with 10 markers
map <- sim.map(c(100,90,80,40), 10, include.x=FALSE)
plotMap(map)

# equally spaced markers
map2 <- sim.map(c(100,90,80,40), 10, include.x=FALSE, eq.spacing=TRUE)
plot(map2)
}

\seealso{ \code{\link{sim.cross}}, \code{\link{plotMap}},
  \code{\link{replace.map}}, \code{\link{pull.map}} }

\keyword{datagen}
