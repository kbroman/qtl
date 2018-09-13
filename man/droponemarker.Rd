\name{droponemarker}
\alias{droponemarker}

\title{Drop one marker at a time and determine effect on genetic map}

\description{
  Drop one marker at a time from a genetic map and calculate the change
  in log likelihood and in the chromosome length, in order to identify
  problematic markers.
}

\usage{
droponemarker(cross, chr, error.prob=0.0001,
                map.function=c("haldane","kosambi","c-f","morgan"),
                m=0, p=0, maxit=4000, tol=1e-6, sex.sp=TRUE,
                verbose=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
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
 \item{verbose}{If TRUE, print information on progress; if > 1, print
   even more information.}
}

\value{
  A data frame (actually, an object of class \code{"scanone"}, so that
  one may use \code{\link{plot.scanone}},
  \code{\link{summary.scanone}}, etc.) with each row being a marker.
  The first two columns are the chromosome ID and position.  The third
  column is a LOD score comparing the hypothesis that the marker is not
  linked to the hypothesis that it belongs at that position.

  In the case of a 4-way cross, with \code{sex.sp=TRUE}, there are two
  additional columns with the change in the estimated female and male genetic lengths
  of the respective chromosome, upon deleting that marker.  
  With \code{sex.sp=FALSE}, or for other types of crosses, there is one
  additional column, with the change in estimated genetic length of the respective
  chromosome, when the marker is omitted.

  A well behaved marker will have a negative LOD score and a small
  change in estimated genetic length.  A poorly behaved marker will have a large
  positive LOD score and a large change in estimated genetic length.  But note
  that dropping the first or last marker on a chromosome could result in
  a large change in estimated length, even if they are not badly
  behaved; for these markers one should focus on the LOD scores, with a large
  positive LOD score being bad.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.bc)
droponemarker(fake.bc, 7, error.prob=0, verbose=FALSE)
}


\seealso{ \code{tryallpositions}, \code{\link{est.map}}, \code{\link{ripple}},
  \code{\link{est.rf}}, \code{\link{switch.order}},
  \code{\link{movemarker}}, \code{\link{drop.markers}} }

\keyword{utilities}
