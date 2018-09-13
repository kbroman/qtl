\name{convert2risib}
\alias{convert2risib}

\title{Convert a cross to RIL by sib mating}

\description{
  Convert a cross to type \code{"risib"} (RIL by sib mating).
}

\usage{
convert2risib(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
}

\value{
  The input cross object, with genotype codes possibly changed and cross
  type changed to \code{"risib"}.
}

\details{
  If there are more genotypes with code 3 (BB) than code 2 (AB), we omit
  the genotypes with code==2 and call those with code==3 the BB
  genotypes.

  If, instead, there are more genotypes with code 2 than code 3, we omit
  the genotypes with code==3 and call those with code==2 the BB genotypes.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(hyper)
hyper.as.risib <- convert2risib(hyper)
}

\seealso{ \code{\link{convert2riself}} }

\keyword{manip}
