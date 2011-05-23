\name{typingGap}
\alias{typingGap}

\title{Maximum distance between genotyped markers}

\description{
  Calculates, for each individual on each chromosome, the maximum
  distance between genotyped markers.
}

\usage{
typingGap(cross, chr)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
   This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.}
}

\details{
  We consider not just the distances between internal genotypes, but
  also distances from the beginning of the chromosome to the first typed
  marker, and similarly for the end of the chromosome.  (The start and end
  of a chromosome are taken to be the locations of the initial and final
  markers.)  
}

\value{
  A matrix with rows corresponding to individuals and columns
  corresponding to chromosomes.  (If there is just one chromosome, it is
  a numeric vector rather than a matrix.)
}

\examples{
data(hyper)
plot(typingGap(hyper, chr=5), ylab="Maximum gap between typed markers (cM)",
     ylim=c(0, diff(range(pull.map(hyper,chr=5)[[1]]))))
}

\seealso{ \code{\link{ntyped}},
          \code{\link{nmissing}}, \code{\link{locateXO}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{utilities}
