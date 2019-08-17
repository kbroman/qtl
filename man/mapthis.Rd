\name{mapthis}
\alias{mapthis}
\docType{data}

\title{Simulated data for illustrating genetic map construction}

\description{
  Simulated data for an F2 intercross, obtained using
  \code{\link{sim.cross}}, useful for illustrating the process of
  constructing a genetic map.
}

\usage{data(mapthis)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details.
}

\details{ These are simulated data, consisting of 300 F2 individuals
typed at 100 markers on five chromosomes. There are no real phenotypes,
just a set of individual identifiers.  The data were simulated for the
purpose of illustrating the process of constructing a genetic map.  The
markers are all assigned to a single chromosome and in a random order,
and there are a number of problematic markers and individuals.

See \url{https://rqtl.org/tutorials/geneticmaps.pdf} for a tutorial
on how to construct a genetic map with these data.}


\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(mapthis)
summary(mapthis)
plot(mapthis)
}

\references{
Broman, K. W. (2010) Genetic map construction with R/qtl.  Technical
report #214, Department of Biostatistics and Medical Informatics,
University of Wisconsin--Madison
}


\seealso{ \code{\link{fake.f2}}, \code{\link{est.rf}},
  \code{\link{est.map}}, \code{\link{formLinkageGroups}},
  \code{\link{orderMarkers}} }

\keyword{datasets}
