\name{convert2riself}
\alias{convert2riself}

\title{Convert a cross to RIL by selfing}

\description{
  Convert a cross to type \code{"riself"} (RIL by selfing).
}

\usage{
convert2riself(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
}

\value{
  The input cross object, with genotype codes possibly changed and cross
  type changed to \code{"riself"}.
}

\details{
  If there are more genotypes with code 3 (BB) than code 2 (AB), we omit
  the genotypes with code==2 and call those with code==3 the BB
  genotypes.

  If, instead, there are more genotypes with code 2 than code 3, we omit
  the genotypes with code==3 and call those with code==2 the BB genotypes.

  Any chromosomes with class \code{"X"} (X chromosome) are changed to
  class \code{"A"} (autosomal).
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
hyper.as.riself <- convert2riself(hyper)
}

\seealso{ \code{\link{convert2risib}} }

\keyword{manip}
