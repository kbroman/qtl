\name{tryallpositions}
\alias{tryallpositions}

\title{Test all possible positions for a marker}

\description{
  Try all possible positions for a marker, keeping all other markers
  fixed, and evaluate the log likelihood and estimate the chromosome
  length. 
}

\usage{
tryallpositions(cross, marker, chr, error.prob=0.0001,
                map.function=c("haldane","kosambi","c-f","morgan"),
                m=0, p=0, maxit=4000, tol=1e-6, sex.sp=TRUE,
                verbose=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{marker}{Character string with name of the marker to move about.}
 \item{chr}{A vector specifying which chromosomes to test for the
   position of the marker.    This should be a vector of character
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
  one may use \code{\link[qtl]{plot.scanone}},
  \code{\link[qtl]{summary.scanone}}, etc.) with each row being a
  possible position for the marker. 
  The first two columns are the chromosome ID and position.  The third
  column is a LOD score comparing the hypotheses that the marker in that
  position versus the hypothesis that it is not linked to that chromosome.

  In the case of a 4-way cross, with \code{sex.sp=TRUE}, there are two
  additional columns with the estimated female and male genetic lengths
  of the respective chromosome, when the marker is in that position.
  With \code{sex.sp=FALSE}, or for other types of crosses, there is one
  additional column, with the estimated genetic length of the respective
  chromosome, when the marker is in that position.

  The row names indicate the nearest flanking markers for each interval.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.bc)
tryallpositions(fake.bc, "D7M301", 7, error.prob=0, verbose=FALSE)
}


\seealso{ \code{\link[qtl]{est.map}}, \code{\link[qtl]{ripple}},
  \code{\link[qtl]{est.rf}}, \code{\link[qtl]{switch.order}},
  \code{\link[qtl]{movemarker}} }

\keyword{utilities}
