\name{jittermap}
\alias{jittermap}

\title{Jitter marker positions in a genetic map}

\description{
  Jitter the marker positions in a genetic map so that no two markers
  are on top of each other.
}

\usage{jittermap(object, amount=1e-6)}

\arguments{
 \item{object}{Either a cross (an object of class \code{cross}; see
   \code{\link[qtl]{read.cross}} for details) or a map (an object of
   class \code{map}; see \code{\link[qtl]{pull.map}} for details).}
 \item{amount}{The amount by which markers should be moved.}
}

\value{
  Either the input cross object or the input map, but with marker
  positions slightly jittered.  If the input was a cross, the function
  \code{\link[qtl]{clean}} is run to strip off any intermediate
  calculations. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(hyper)
hyper <- jittermap(hyper)
}

\seealso{ \code{\link[qtl]{pull.map}}, \code{\link[qtl]{replace.map}},
  \code{\link[qtl]{summary.cross}} }

\keyword{manip}
