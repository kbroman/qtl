\name{scanonevar.meanperm}
\alias{scanonevar.meanperm}

\title{Permutation test for mean effect in scanonevar}

\description{Executes permutations of the genotypes in the mean-effect part of scanonevar}

\usage{
scanonevar.meanperm(cross, pheno.col=1, mean_covar=NULL, var_covar=NULL,
           n.mean.perm = 2, seed = 27517, maxit=25, tol=1e-6)
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
  \item{n.mean.perm}{Numeric vector of length one indicates the number of permutations to execute.}
  \item{seed}{Numeric vector of length one indicates the random seed to start the permutations.}
  \item{maxit}{Maximum number of iterations in the algorithm to fit the
    model at a given position.}
  \item{tol}{Tolerance for convergence.}
  \item{quiet}{If \code{FALSE}, print some information about the course
    of the calculations.}
}


\value{
  A vector of length \code{n.mean.perm} of the maximum negative log10 p-value that resulted from each permutation.
}
