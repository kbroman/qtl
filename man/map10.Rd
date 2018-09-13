\name{map10}
\alias{map10}
\docType{data}

\title{An example genetic map}

\description{
  A genetic map corresponding approximately to the mouse genome with
  a 10 cM marker spacing.
}

\usage{data(map10)}

\format{
  An object of class \code{map}: a list whose components are vectors of
  marker locations.  This map approximates the mouse genome, with 20
  chromosomes (including the X chromosome) and 187 markers at an
  approximately 10 cM spacing.  The markers are equally spaced on each
  chromosome, but the spacings are a bit above or below 10 cM, so that
  the lengths match those in the Mouse Genome Database.}

\examples{
  data(map10)
  plot(map10)

  mycross <- sim.cross(map10, type="f2", n.ind=100)
}

\seealso{ \code{\link{sim.map}}, \code{\link{plotMap}},
  \code{\link{pull.map}} }

\keyword{datasets}
