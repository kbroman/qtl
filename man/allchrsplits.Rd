\name{allchrsplits}
\alias{allchrsplits}

\title{Test all possible splits of a chromosome into two pieces}

\description{
  In order to assess the support for a linkage group, this function
  splits the linkage groups into two pieces at each interval and in each
  case calculates a LOD score comparing the combined linkage group to
  the two pieces.
}

\usage{
allchrsplits(cross, chr, error.prob=0.0001,
                map.function=c("haldane","kosambi","c-f","morgan"),
                m=0, p=0, maxit=4000, tol=1e-6, sex.sp=TRUE,
                verbose=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{A vector specifying which chromosomes to study.
   This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.} 
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions. (Ignored if m > 0.)}
 \item{m}{Interference parameter for the chi-square model for
   interference; a non-negative integer, with m=0 corresponding to no
   interference. This may be used only for a backcross or intercross.}
 \item{p}{Proportion of chiasmata from the NI mechanism, in the Stahl
   model; p=0 gives a pure chi-square model.  This may be used only for
   a backcross or intercross.} 
 \item{maxit}{Maximum number of EM iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps; this is 
 used only for the 4-way cross.}
 \item{verbose}{If TRUE, print information on progress.}
}

\value{
  A data frame (actually, an object of class \code{"scanone"}, so that
  one may use \code{\link{plot.scanone}},
  \code{\link{summary.scanone}}, etc.) with each row being an interval
  at which a split is made.
  The first two columns are the chromosome ID and midpoint of the interval.  The third
  column is a LOD score comparing the combined linkage group to the
  split into two linkage groups.  A fourth column (\code{gap}) indicates the length of
  each interval.

  The row names indicate the flanking markers for each interval.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
allchrsplits(fake.bc, 7, error.prob=0, verbose=FALSE)
}


\seealso{ \code{\link{est.map}}, \code{\link{ripple}},
  \code{\link{est.rf}}, \code{\link{switch.order}},
  \code{\link{movemarker}} }

\keyword{utilities}
