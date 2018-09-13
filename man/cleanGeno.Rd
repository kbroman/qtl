\name{cleanGeno}
\alias{cleanGeno}

\title{Delete genotypes that are possibly in error}

\description{
  Delete genotypes from a cross that are indicated to be possibly in
  error, as they result in apparent tight double-crossovers.
}

\usage{
cleanGeno(cross, chr, maxdist=2.5, maxmark=2, verbose=TRUE)
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
 \item{maxdist}{A vector specifying the maximum distance between two
 crossovers.} 
 \item{maxmark}{A vector specifying the maximum number of typed markers
 between two crossovers.}
 \item{verbose}{If TRUE, print information on the numbers of genotypes
 omitted from each chromosome.}

}

\details{
  We first use \code{\link{locateXO}} to identify crossover locations.
  If a pair of adjacted crossovers are separated by no more than
  \code{maxdist} and contain no more than \code{maxmark} genotyped
  markers, the intervening genotypes are omitted (that is, changed to
  \code{NA}).  

  The arguments \code{maxdist} and \code{maxmark} may be vectors.  (If
  both have length greater than 1, they must have the same length.)  If
  they are vectors, genotypes are omitted if they satisify any one of
  the (\code{maxdist}, \code{maxmark}) pairs.
}

\value{
  The input \code{cross} object with suspect genotypes omitted.
}

\examples{
data(hyper)
sum(ntyped(hyper))
hyperc <- cleanGeno(hyper, chr=4, maxdist=c(2.5, 10), maxmark=c(2, 1))
sum(ntyped(hyperc))
}

\seealso{ \code{\link{locateXO}},
          \code{\link{countXO}}, \code{\link{calc.errorlod}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{utilities}
