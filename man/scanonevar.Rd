\name{scanonevar}
\alias{scanonevar}

\title{Genome scan for QTL affecting mean and/or variance}

\description{
  Genome scan with a single QTL model using a double generalized linear
  model, for loci that can affect the variance as well as the mean.
}

\usage{
scanonevar(cross, pheno.col=1, mean_covar=NULL, var_covar=NULL,
           family=gaussian(), maxit=25,
           fixed_gamma_disp=FALSE, quiet=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This must be a single value (integer index or
   phenotype name) or a numeric vector of
   phenotype values, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
  \item{mean_covar}{Numeric matrix with covariates affecting the mean.}
  \item{var_covar}{Numeric matrix with covariates affecting the variances.}
  \item{family}{Family object for \code{\link[dglm]{dglm}}.}
  \item{maxit}{Maximum number of iterations in the algorithm to fit the
    model at a given position.}
  \item{fixed_gamma_disp}{If \code{FALSE}, estimate the dispersion
    parameter in the gamma GLM for the residual variance.}
  \item{quiet}{If \code{FALSE}, print some information about the course
    of the calculations.}
}


\value{
  A data frame (with class \code{"scanone"}, in the form output by
  \code{\link{scanone}}), with four columns: chromosome, position, the -log P-value for
  the mean effect, and the -log P-value for the effect on the variance.
  The result is given class \code{"scanone"}
}

\references{
  Dunn, P. K. and Smyth, G. K. (2012) dglm: Double generalized linear
  models, \url{http://CRAN.R-project.org/package=dglm}.

  Ronnegard, L. and Valdar W. (2011) Detecting major genetic loci
  controlling phenotypic variability in experimental crosses. Genetics
  188:435-447
}

\author{Lars Ronnegard and Karl Broman}

\examples{
data(fake.bc)
\dontshow{fake.bc <- fake.bc[18:19,]}
fake.bc <- calc.genoprob(fake.bc, step=2.5)
out <- scanonevar(fake.bc)
color <- c("slateblue", "violetred")
plot(out, lod=1:2, col=color, bandcol="gray80")
legend("topright", lwd=2, c("mean", "variance"), col=color)

# use format="allpeaks" to get summary for each of mean and variance
#  also consider format="tabByCol" or format="tabByChr"
summary(out, format="allpeaks")

# with sex and age as covariates
covar <- fake.bc$pheno[,c("sex", "age")]
out.cov <- scanonevar(fake.bc, mean_covar=covar, var_covar=covar)
}

\seealso{ \code{\link[dglm]{dglm}}, \code{\link{scanone}},
  \code{\link{summary.scanone}}, \code{\link{calc.genoprob}},
  \code{\link{summary.scanoneperm}} }

\keyword{models}
