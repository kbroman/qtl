\name{est.rf}
\alias{est.rf}

\title{Estimate pairwise recombination fractions}

\description{
  Estimate the sex-averaged recombination fraction between all pairs of
  genetic markers.
}

\usage{
est.rf(cross, maxit=10000, tol=1e-6)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{maxit}{Maximum number of iterations for the EM algorithm (not
   used with backcrosses).}
 \item{tol}{Tolerance for determining convergence (not used with
   backcrosses).} 
}

\value{
  The input \code{cross} object is returned with a component, \code{rf},
  added.  This is a matrix of size (tot.mar x tot.mar).  The diagonal
  contains the number of typed meioses per marker, the lower triangle
  contains the estimated recombination fractions, and the upper triangle
  contains the LOD scores (testing rf = 0.5).
}

\details{
  For a backcross, one can simply count recombination events.  For an
  intercross or 4-way cross, a version of the EM algorithm must be used
  to estimate recombination fractions. (Since, for example, in an
  intercross individual that is heterozygous at two loci, it is not
  known whether there were 0 or 2 recombination events.)  Note that, for
  the 4-way cross, we estimate sex-averaged recombination fractions.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(badorder)
badorder <- est.rf(badorder)
plotRF(badorder)
}

\seealso{ \code{\link{plotRF}}, \code{\link{pull.rf}}, \code{\link{plot.rfmatrix}},
  \code{\link{est.map}},
  \code{\link{badorder}}, \code{\link{checkAlleles}} }

\keyword{utilities}
