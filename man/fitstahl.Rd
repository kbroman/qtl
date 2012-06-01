\name{fitstahl}
\alias{fitstahl}

\title{Fit Stahl interference model}

\description{
  Fit the Stahl model for crossover inference (or the chi-square model,
  which is a special case).
}

\usage{
fitstahl(cross, chr, m, p, error.prob=0.0001, maxit=4000, tol=1e-4,
         maxm=15, verbose=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
\item{m}{Interference parameter (a non-negative integer); if
   unspecified, this is estimated.}
 \item{p}{The proportion of chiasmata coming from the no interference
   mechanism in the Stahl model (0 <= p <= 1).  p=0 gives the chi-square
   model.  If unspecified, this is estimated.}
 \item{error.prob}{The genotyping error probability.  If = NULL, it is
   estimated.}
 \item{maxit}{Maximum number of iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{maxm}{Maximum value of m to consider, if m is unspecified.}
 \item{verbose}{Logical; indicates whether to print tracing information.}
}

\details{
  This function is currently only available for backcrosses and
  intercrosses. 

  The Stahl model of crossover interference (of which the chi-square
  model is a special case) is fit.      In the chi-square model, points
  are tossed down onto the four-strand bundle according to a Poisson
  process, and every \eqn{(m+1)}st point is a chiasma.  With the
  assumption of no chromatid interference, crossover locations on a
  random meiotic product are obtained by thinning the chiasma process.
  The parameter \eqn{m} (a non-negative integer) governs the strength of
  crossover interference, with \eqn{m=0} corresponding to no
  interference. 

  In the Stahl model, chiasmata on the four-strand bundle are a
  superposition of chiasmata from two mechanisms, one following a
  chi-square model and one exhibiting no interference.  An additional
  parameter, \eqn{p}, gives the proportion of chiasmata from the no
  interference mechanism.

  If all of \code{m}, \code{p}, and \code{error.prob} are specified, any
  of them with length > 1 must all have the same length.  

  If \code{m} is unspecified, we do a grid search starting at 0 and stop
  when the likelihood decreases (thus assuming a single mode), or
  \code{maxm} is reached.
}

\value{
  A matrix with four columns: m, p, error.prob, and the log likelihood.

  If specific values for m, p, error.prob are provided, the log
  likelihood for each set are given.

  If some are left unspecified, the maximum likelihood estimates are
  provided in the results.  
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
\dontshow{set.seed(12201969)}
# Simulate genetic map: one chromosome of length 200 cM with
# a 2 cM marker spacing
mymap <- sim.map(200, 51, anchor.tel=TRUE, include.x=FALSE,
                 sex.sp=FALSE, eq.spacing=TRUE)

# Simulate data under the chi-square model, no errors
mydata <- sim.cross(mymap, n.ind=250, type="bc",
                    error.prob=0, m=3, p=0)

# Fit the chi-square model for specified m's
\dontrun{output <- fitstahl(mydata, m=1:5, p=0, error.prob=0)
}\dontshow{output <- fitstahl(mydata, m=1, p=0, error.prob=0,tol=0.001)}
plot(output$m, output$loglik, lwd=2, type="b")

# Find the MLE of m in the chi-square model
\dontrun{mle <- fitstahl(mydata, p=0, error.prob=0)}

\dontrun{
# Simulate data under the Stahl model, no errors
mydata <- sim.cross(mymap, n.ind=250, type="bc",
                    error.prob=0, m=3, p=0.1)

# Find MLE of m for the Stahl model with known p 
mle.stahl <- fitstahl(mydata, p=0.1, error.prob=0)

# Fit the Stahl model with unknown p and m,
# get results for m=0, 1, 2, ..., 8
output <- fitstahl(mydata, m=0:8, error.prob=0)
plot(output$m, output$loglik, type="b", lwd=2)}
}

\references{
  Armstrong, N. J., McPeek, M. J. and Speed, T. P. (2006) Incorporating
  interference into linkage analysis for experimental crosses.
  \emph{Biostatistics} \bold{7}, 374--386.
  
  Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis of
  crossover interference using the chi-square model.  \emph{Genetics}
  \bold{139}, 1045--1056.
}  

\seealso{ \code{\link{est.map}}, \code{\link{sim.cross}}}

\keyword{models}
