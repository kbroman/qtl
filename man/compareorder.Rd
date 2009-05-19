\name{compareorder}
\alias{compareorder}

\title{Compare two orderings of markers on a chromosome}

\description{
  Compare the likelihood of an alternative order for markers on a
  chromosome to the current order.
}

\usage{
compareorder(cross, chr, order, error.prob=0.0001, 
             map.function=c("haldane","kosambi","c-f","morgan"),
             maxit=4000, tol=1e-6, sex.sp=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
\item{chr}{The chromosome to investigate.  Only one chromosome is
  allowed. (This should be a character string referring to the
  chromosomes by name.)} 
\item{order}{The alternate order of markers on the chromosome: a numeric
             vector that is a permutation of the integers from 1 to the
             number of markers on the chromosome.}
 \item{error.prob}{Assumed genotyping error rate used in the calculation
   of the penetrance Pr(observed genotype | true genotype).}
 \item{map.function}{Indicates whether to use the Haldane, Kosambi,
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
 \item{maxit}{Maximum number of EM iterations to perform.}
 \item{tol}{Tolerance for determining convergence.}
 \item{sex.sp}{Indicates whether to estimate sex-specific maps; this is 
   used only for the 4-way cross.}
}

\value{
  A data frame with two rows: the current order in the input cross object, and the
  revised order.  The first column is the log10 likelihood of the new order relative
  to the original one (positive values indicate that the new order is better supported). 
  The second column is the estimated genetic length of the chromosome for each order.  
  In the case of sex-specific maps, there are separate columns for the female and male
  genetic lengths. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(badorder)
compareorder(badorder, chr=1, order=c(1:8,11,10,9,12))
}

\seealso{ \code{\link[qtl]{ripple}}, \code{\link[qtl]{switch.order}},
  \code{\link[qtl]{move.marker}} }

\keyword{utilities}
