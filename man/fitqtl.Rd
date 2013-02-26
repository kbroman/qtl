\name{fitqtl}
\alias{fitqtl}

\title{Fit a multiple-QTL model}

\description{
  Fits a user-specified multiple-QTL model.  If specified,
  a drop-one-term analysis will be performed.
}

\usage{
fitqtl(cross, pheno.col=1, qtl, covar=NULL, formula, method=c("imp", "hk"),
       model=c("normal", "binary"), dropone=TRUE, get.ests=FALSE,
       run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
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
   in the data frame \code{covar}.  }
 \item{method}{Indicates whether to use multiple imputation or
   Haley-Knott regression.}
 \item{model}{The phenotype model: the usual model or a model for binary
  traits}
 \item{dropone}{If TRUE, do drop-one-term analysis.}
 \item{get.ests}{If TRUE, return estimated QTL effects and their
   estimated variance-covariance matrix.}
 \item{run.checks}{If TRUE, check the input formula and check for
   individuals with missing phenotypes or covariates.}
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

  In the drop-one-term analysis, for a given QTL/covariate model, all
  submodels will be analyzed.  For each term in the input formula, when
  it is dropped, all higher order terms that contain it will also be
  dropped.  The comparison between the new model and the full (input)
  model will be output.

  The estimated percent variances explained for the QTL are simplify
  transformations of the conditional LOD scores by the formula \eqn{h^2
  = 1 - 10^{-(2/n) {\rm LOD}}}.  While these may be reasonable for
  unlinked, additive QTL, \bold{they can be completely wrong in the case
  of linked QTL}, but we don't currently have any alternative.

  For \code{model="binary"}, a logistic regression model is used.

  \bold{The part to get estimated QTL effects is fully working only for
    the case of autosomes in a backcross, intercross, RIL or doubled haploids.
    In other cases the values returned are based on a
    design matrix that is convenient for calculations but not easily
    interpreted.}

  The estimated QTL effects for a backcross are derived by the coding
  scheme \eqn{\pm}{+/-} 1/2 for AA and AB, so that the additive
  effect corresponds to the difference between phenotype averages for
  the two genotypes.  For doubled haploids and RIL, the coding scheme is
  \eqn{\pm}{+/-} 1 for AA and BB, so that the additive effect
  corresponds to half the difference between the phenotype averages for
  the two homozygotes.

  For an intercross, the additive effect is derived from the coding
  scheme -1/0/+1 for genotypes AA/AB/BB, and so is half the difference
  between the phenotype averages for the two homozygotes.  The dominance
  deviation is derived from the coding scheme 0/+1/0 for genotypes
  AA/AB/BB, and so is the difference between the phenotype average for
  the heterozygotes and the midpoint between the phenotype averages for
  the two homozygotes.

  Epistatic effects and QTL \eqn{\times}{x} covariate interaction
  effects are obtained through the products of the corresponding
  additive/dominant effect columns.
}

\value{
  An object of class \code{fitqtl}. It may contains as many as four components:
  \itemize{
    \item \code{result.full} is the ANOVA table as a matrix for the full model
    result. It contains the degree of freedom (df), Sum of squares (SS),
    mean square (MS), LOD score (LOD), percentage of variance explained
    (\%var) and P value (Pvalue).

    \item \code{lod} is the LOD score from the fit of the full model.

    \item \code{result.drop} is a drop-one-term ANOVA table as a
    matrix. It contains degrees of freedom (df), Type III sum of squares
    (Type III SS), LOD score(LOD), percentage of variance explained
    (\%var), F statistics (F value), and P values for chi square
    (Pvalue(chi2)) and F distribution (Pvalue(F)).  Note that the degree
    of freedom, Type III sum of squares, the LOD score and the
    percentage of variance explained are the values comparing the full
    to the sub-model with the term dropped. Also note that for
    imputation method, the percentage of variance explained, the the F
    values and the P values are approximations calculated from the LOD
    score.

    \item \code{ests} contains the estimated QTL effects and standard errors.
}

  \bold{The part to get estimated QTL effects is fully working only for
    the case of autosomes in a backcross, intercross, RIL or doubled haploids.
    In other cases the values returned are based on a
    design matrix that is convenient for calculations but not easily
    interpreted.}
}


\references{
  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Hao Wu; Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.f2)

# take out several QTLs and make QTL object
qc <- c(1, 8, 13)
qp <- c(26, 56, 28)
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

# fit model with 3 interacting QTLs interacting
# (performing a drop-one-term analysis)
lod <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1*Q2*Q3, method="hk")
summary(lod)

\dontrun{
# fit an additive QTL model
lod.add <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1+Q2+Q3, method="hk")
summary(lod.add)

# fit the model including sex as an interacting covariate
Sex <- data.frame(Sex=pull.pheno(fake.f2, "sex"))
lod.sex <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1*Q2*Q3*Sex,
                  cov=Sex, method="hk")
summary(lod.sex)

# fit the same with an additive model
lod.sex.add <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1+Q2+Q3+Sex,
                      cov=Sex, method="hk")
summary(lod.sex.add)
}}

\seealso{ \code{\link{summary.fitqtl}}, \code{\link{makeqtl}},
\code{\link{scanqtl}}, \code{\link{refineqtl}},
  \code{\link{addtoqtl}} ,
  \code{\link{dropfromqtl}},
  \code{\link{replaceqtl}},
  \code{\link{reorderqtl}}}

\keyword{models}
