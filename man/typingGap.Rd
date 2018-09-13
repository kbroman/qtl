\name{typingGap}
\alias{typingGap}

\title{Maximum distance between genotyped markers}

\description{
  Calculates, for each individual on each chromosome, the maximum
  distance between genotyped markers.
}

\usage{
typingGap(cross, chr, terminal=FALSE)
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
\item{terminal}{If TRUE, just look at terminal typing gaps (from the
 terminal markers to the first typed marker).}
}

\details{
  We consider not just the distances between internal genotypes, but
  also distances from the beginning of the chromosome to the first typed
  marker, and similarly for the end of the chromosome.  (The start and end
  of a chromosome are taken to be the locations of the initial and final
  markers.)  If \code{terminal=TRUE}, we look only at those beginning
  and end distances.
}

\value{
  A matrix with rows corresponding to individuals and columns
  corresponding to chromosomes.  (If there is just one chromosome, it is
  a numeric vector rather than a matrix.)
}

\examples{
data(hyper)
plot(typingGap(hyper, chr=5),
     ylab="Maximum gap between typed markers (cM)",
     ylim=c(0, diff(range(pull.map(hyper,chr=5)[[1]]))))

plot(typingGap(hyper, chr=4),
     ylab="Maximum gap between typed markers (cM)",
     ylim=c(0, diff(range(pull.map(hyper,chr=4)[[1]]))))

plot(typingGap(hyper, chr=4, terminal=TRUE),
     ylab="Maximum gap between chr end and typed marker (cM)",
     ylim=c(0, diff(range(pull.map(hyper,chr=4)[[1]]))))
}

\seealso{ \code{\link{ntyped}},
          \code{\link{nmissing}}, \code{\link{locateXO}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{utilities}
