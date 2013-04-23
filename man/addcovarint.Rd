\name{addcovarint}
\alias{addcovarint}

\title{Add QTL x covariate interaction to a multiple-QTL model}

\description{
  Try adding all QTL x covariate interactions, one at a time, to a
  multiple QTL model, for a given set of covariates.
}

\usage{
addcovarint(cross, pheno.col=1, qtl, covar=NULL, icovar, formula, 
            method=c("imp","hk"), model=c("normal", "binary"), 
            verbose=TRUE, pvalues=TRUE, simple=FALSE, tol=1e-4,
            maxit=1000, require.fullrank=TRUE)
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
 \item{icovar}{Vector of character strings indicating the columns in
   \code{covar} to be considered for QTL x covariate interactions.}
 \item{formula}{An object of class \code{\link[stats]{formula}}
   indicating the model to be fitted.  (It can also be the character
   string representation of a formula.)  QTLs are referred to as
   \code{Q1}, \code{Q2}, etc.  Covariates are referred to by their names
   in the data frame \code{covar}.} 
 \item{method}{Indicates whether to use multiple imputation or
   Haley-Knott regression.}
 \item{model}{The phenotype model: the usual model or a model for binary
  traits} 
 \item{verbose}{If TRUE, will print a message if there are no
   interactions to test.}
 \item{pvalues}{If FALSE, p-values will not be included in the results.}
 \item{simple}{If TRUE, don't include p-values or sums of squares in the summary.}
 \item{tol}{Tolerance for convergence for the binary trait model.}
 \item{maxit}{Maximum number of iterations for fitting the binary trait
   model.} 
 \item{require.fullrank}{If TRUE, give LOD=0 when covariate matrix in
   the linear regression is not of full rank.}
}

\details{
  The formula is used to specified the model to be fit. In the
  formula, use \code{Q1}, \code{Q2}, etc., or \code{q1},
  \code{q2}, etc., to represent the QTLs, and the column names in the
  covariate data frame to represent the covariates.
    
  We enforce a hierarchical structure on the model formula: if a QTL or
  covariate is in involved in an interaction, its main effect must also
  be included.   
}

\section{Value}{
  An object of class \code{addcovarint}, with results as in the
  drop-one-term analysis from \code{\link{fitqtl}}.  This is a data
  frame (given class \code{"addcovarint"}, with the following columns: 
  degrees of freedom (df), Type III sum of squares (Type III
  SS), LOD score(LOD), percentage of variance explained (\%var), F
  statistics (F value),  and P values for chi square (Pvalue(chi2))
  and F distribution (Pvalue(F)).

  Note that the degree of freedom, Type III sum of squares, the LOD
  score and the percentage of variance explained are the values
  comparing the full to the sub-model with the term dropped. Also note
  that for imputation method, the percentage of variance explained, the
  the F values and the P values are approximations calculated from the
  LOD score.

  QTL x covariate interactions already included in the input \code{formula} are
  not tested. 
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
fake.f2 <- subset(fake.f2, chr=qc)
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

# use the sex phenotype as the covariate
covar <- data.frame(sex=fake.f2$pheno$sex)

# try all possible QTL x sex interactions, one at a time
addcovarint(fake.f2, pheno.col=1, qtl, covar, "sex", y~Q1+Q2+Q3,
            method="hk")
}

\seealso{ \code{\link{addint}}, \code{\link{fitqtl}}, \code{\link{makeqtl}},
  \code{\link{scanqtl}}, \code{\link{refineqtl}},
  \code{\link{addqtl}}, \code{\link{addpair}} }

\keyword{models}
