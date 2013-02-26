\name{addqtl}
\alias{addqtl}

\title{Scan for an additional QTL in a multiple-QTL model}

\description{
  Scan for an additional QTL in the context of a
  multiple QTL model.
}

\usage{
addqtl(cross, chr, pheno.col=1, qtl, covar=NULL, formula,
       method=c("imp","hk"), model=c("normal", "binary"),
       incl.markers=TRUE, verbose=FALSE, tol=1e-4, maxit=1000,
       forceXcovar=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to be scanned. If
   missing, all chromosomes are scanned. Refer to chromosomes by
   name. Refer to chromosomes with a preceding \code{-} to have all
   chromosomes but those considered.  A logical (TRUE/FALSE) vector may
   also be used.}
 \item{pheno.col}{Column number in the phenotype matrix to be
   used as the phenotype.  One may also give a character string matching
   a phenotype name. Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{qtl}{An object of class \code{qtl}, as output from
   \code{\link{makeqtl}}.}
 \item{covar}{A matrix or data.frame of covariates.  These must be
   strictly numeric.}
 \item{formula}{An object of class \code{\link[stats]{formula}}
   indicating the model to be fitted.  (It can also be the character
   string representation of a formula.)  QTLs are referred to as
   \code{Q1}, \code{Q2}, etc.  Covariates are referred to by their names
   in the data frame \code{covar}.  If the new QTL is not included in
   the formula, its main effect is added.}
 \item{method}{Indicates whether to use multiple imputation or
   Haley-Knott regression.}
 \item{model}{The phenotype model: the usual model or a model for binary
  traits}
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.  If \code{\link{calc.genoprob}} or
   \code{\link{sim.geno}} were run with
   \code{stepwidth="variable"} or \code{stepwidth="max"}, we force \code{incl.markers=TRUE}.}
 \item{verbose}{If TRUE, display information about the progress of
   calculations.  If \code{verbose} is an integer > 1, further messages
   from \code{\link{scanqtl}} are also displayed.}
 \item{tol}{Tolerance for convergence for the binary trait model.}
 \item{maxit}{Maximum number of iterations for fitting the binary trait
   model.}
 \item{forceXcovar}{If TRUE, force inclusion of X-chr-related covariates
   (like sex and cross direction).}
}

\details{
  The formula is used to specified the model to be fit. In the
  formula, use \code{Q1}, \code{Q2}, etc., or \code{q1},
  \code{q2}, etc., to represent the QTLs, and the column names in the
  covariate data frame to represent the covariates.

  We enforce a hierarchical structure on the model formula: if a QTL or
  covariate is in involved in an interaction, its main effect must also
  be included.

  If one wishes to scan for QTL that interact with another QTL, include
  it in the formula (with an index of one more than the number of QTL in
  the input \code{qtl} object).
}

\section{Value}{
  An object of class \code{scanone}, as produced by the
  \code{\link{scanone}} function.  LOD scores are relative to the
  base model (with any terms that include the new QTL omitted).
}

\references{
  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 8, 13)
qp <- c(26, 56, 28)

fake.f2 <- subset(fake.f2, chr=c(1,2,3,8,13))
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}

fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

# scan for an additional QTL
out1 <- addqtl(fake.f2, qtl=qtl, formula=y~Q1+Q2+Q3, method="hk")
max(out1)

# scan for an additional QTL that interacts with the locus on chr 1
out2 <- addqtl(fake.f2, qtl=qtl, formula=y~Q1*Q4+Q2+Q3, method="hk")
max(out2)

# plot interaction LOD scores
plot(out2-out1)
}

\seealso{ \code{\link{scanone}}, \code{\link{fitqtl}},
  \code{\link{scanqtl}}, \code{\link{refineqtl}},
  \code{\link{makeqtl}}, \code{\link{addtoqtl}},
  \code{\link{addpair}}, \code{\link{addint}} }

\keyword{models}
