\name{nmissing}
\alias{nmissing}

\title{Number of missing genotypes}

\description{
  Count the number of missing genotypes for each individual or each marker
  in a cross.
}

\usage{nmissing(cross, what=c("ind","mar"))}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{what}{Indicates whether to count missing genotypes for each
   individual or each marker.}
}

\value{
  A vector containing the number of missing genotypes for each
  individual or for each marker.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(listeria)

# plot number of missing genotypes for each individual
plot(nmissing(listeria))

# plot number of missing genotypes for each marker
plot(nmissing(listeria, what="mar"))
}

\seealso{ \code{\link[qtl]{summary.cross}},
          \code{\link[qtl]{nind}}, \code{\link[qtl]{totmar}} }

\keyword{utilities}
