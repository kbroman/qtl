\name{convert2sa}
\alias{convert2sa}

\title{Convert a sex-specific map to a sex-averaged one}

\description{
  Convert a sex-specific map to a sex-averaged one, assuming that the
  female and male maps are actually the same (that is, that the map was
  estimated assuming a common recombination rate in females and males).
}

\usage{
convert2sa(map, tol=1e-4)
}

\arguments{
 \item{map}{A map object with sex-specific locations (but assuming that
 the female and male maps are the same), as output by the function
 \code{\link{est.map}} for a 4-way cross, with argument
 \code{sex.sp=FALSE}.}
  \item{tol}{Tolerance value for inspecting the differences between the
 female and male maps; if they differ by more than this tolerance, a
 warning is issued.}
}

\value{
  A map object, with sex-averaged distances.  
}

\details{
  We pull out just the female marker locations, and give a warning if
  there are large differences between the female and male maps.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.4way)
\dontrun{fake.4way <- subset(fake.4way, chr="-X")}
\dontshow{fake.4way <- subset(fake.4way, chr=18:19)}
nm <- est.map(fake.4way, sex.sp=FALSE)
plot(convert2sa(nm))
}

\seealso{ \code{\link{est.map}}, \code{\link{plotMap}} }

\keyword{manip}
