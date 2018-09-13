\name{scanoneboot}
\alias{scanoneboot}

\title{Bootstrap to get interval estimate of QTL location}

\description{
  Nonparametric bootstrap to get an estimated confidence interval for
  the location of a QTL, in the context of a single-QTL model.
}

\usage{
scanoneboot(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"),
            method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
            addcovar=NULL, intcovar=NULL, weights=NULL,
            use=c("all.obs", "complete.obs"), upper=FALSE,
            ties.random=FALSE, start=NULL, maxit=4000,
            tol=1e-4, n.boot=1000, verbose=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{The chromosome to investigate.  Only one chromosome is
  allowed. (This should be a character string referring to the
  chromosomes by name.)}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  One may also give a character string matching
   a phenotype name.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}  
 \item{model}{The phenotypic model: the usual normal model, a model for
   binary traits, a two-part model or non-parametric analysis}
 \item{method}{Indicates whether to use the EM algorithm, 
   imputation, Haley-Knott regression, the extended Haley-Knott method,
   or marker regression.  Not all methods are available for all models.
   Marker regression is performed either by dropping individuals with
   missing genotypes (\code{"mr"}), or by first filling in missing data
   using a single imputation (\code{"mr-imp"}) or by the Viterbi
   algorithm (\code{"mr-argmax"}).}
 \item{addcovar}{Additive covariates;
   allowed only for the normal and binary models.}
 \item{intcovar}{Interactive covariates (interact with QTL genotype);
   allowed only for the normal and binary models.}
 \item{weights}{Optional weights of individuals.  Should be either NULL
   or a vector of length n.ind containing positive weights.  Used only
   in the case \code{model="normal"}.}
 \item{use}{In the case that multiple phenotypes are selected to be
   scanned, this argument indicates whether to use all individuals, 
   including those missing some phenotypes, or just those individuals
   that have data on all selected phenotypes.}
 \item{upper}{Used only for the two-part model; if true, the
   "undefined" phenotype is the maximum observed phenotype; otherwise,
   it is the smallest observed phenotype.} 
 \item{ties.random}{Used only for the non-parametric "model"; if TRUE,
   ties in the phenotypes are ranked at random.  If FALSE, average ranks
   are used and a corrected LOD score is calculated.}
 \item{start}{Used only for the EM algorithm with the normal model and
   no covariates.  If \code{NULL}, use the usual starting values; if
   length 1, use random initial weights for EM; otherwise, this should
   be a vector of length n+1 (where n is the number of possible
   genotypes for the cross), giving the initial values for EM.}
 \item{maxit}{Maximum number of iterations for methods \code{"em"} and
   \code{"ehk"}.}
 \item{tol}{Tolerance value for determining convergence for methods
   \code{"em"} and \code{"ehk"}.} 
 \item{n.boot}{Number of bootstrap replicates.}
 \item{verbose}{If TRUE, display information about the progress of the
   bootstrap.} 
}

\details{
  We recommend against the use of the bootstrap to derive a confidence
  interval for the location of a QTL; see Manichaikul et al. (2006).
  Use \code{\link{lodint}} or \code{\link{bayesint}} instead.

  The bulk of the arguments are the same as for the
  \code{\link{scanone}} function.  A single chromosome should be
  indicated with the \code{chr} argument; otherwise, we focus on the
  first chromosome in the input \code{cross} object.

  A single-dimensional scan on the relevant chromosome is performed.  We
  further perform a nonparametric bootstrap (sampling individuals \emph{with
  replacement} from the available data, to create a new data set with
  the same size as the input cross; some individuals with be duplicated
  and some omitted).  The same scan is performed with the resampled data;
  for each bootstrap replicate, we store only the location with maximum
  LOD score. 

  Use \code{\link{summary.scanoneboot}} to obtain the desired
  confidence interval.
}

\value{
  A vector of length \code{n.boot}, giving the estimated QTL locations
  in the bootstrap replicates.  The results for the original data are
  included as an attribute, \code{"results"}.  
}

\references{
  Manichaikul, A., Dupuis, J., Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S} and Broman, K. W. (2006) Poor
  performance of bootstrap confidence intervals for the location of a
  quantitative trait locus. \emph{Genetics} \bold{174}, 481--489.
  
  Visscher, P. M., Thompson, R. and Haley, C. S. (1996) Confidence
  intervals in QTL mapping by bootstrap.  \emph{Genetics} \bold{143},
  1013--1020.  
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2, step=1, err=0.001)
\dontrun{bootoutput <- scanoneboot(fake.f2, chr=13, method="hk")}
\dontshow{bootoutput <- scanoneboot(fake.f2, chr=13, method="hk", n.boot=50)}

plot(bootoutput)
summary(bootoutput)
}

\seealso{ \code{\link{scanone}}, \code{\link{summary.scanoneboot}},
  \code{\link{plot.scanoneboot}},
  \code{\link{lodint}}, \code{\link{bayesint}} }

\keyword{utilities}
