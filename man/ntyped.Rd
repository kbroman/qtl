\name{ntyped}
\alias{ntyped}

\title{Number of genotypes}

\description{
  Count the number of genotypes for each individual or each marker
  in a cross.
}

\usage{ntyped(cross, what=c("ind","mar"))}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{what}{Indicates whether to count genotypes for each
   individual or each marker.}
}

\value{
  A vector containing the number of genotypes for each
  individual or for each marker.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)

# plot number of genotypes for each individual
plot(ntyped(listeria))

# plot number of genotypes for each marker
plot(ntyped(listeria, what="mar"))
}

\seealso{ \code{\link{nmissing}},  \code{\link{summary.cross}},
          \code{\link{nind}}, \code{\link{totmar}} }

\keyword{utilities}
