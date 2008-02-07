\name{fitqtl}
\alias{fitqtl}

\title{Fit a multiple-QTL model}

\description{
  Fits a user-specified multiple-QTL model.  If specified, 
  a drop-one-term analysis will be performed.
}

\usage{
fitqtl(cross, pheno.col=1, qtl, covar=NULL, formula, method=c("imp", "hk"),
       dropone=TRUE, get.ests=FALSE, run.checks=TRUE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This can be a vector of integers; for methods
   \code{"hk"} and \code{"imp"} this can be considerably faster than doing
   them one at a time.}
 \item{qtl}{An object of class \code{qtl}, as output from
   \code{\link[qtl]{makeqtl}}.} 
 \item{covar}{A data.frame of covariates.  These must be strictly
   numeric.} 
 \item{formula}{An object of class \code{\link[stats]{formula}}
   indicating the model to be fitted.  QTLs are referred to as
   \code{Q1}, \code{Q2}, etc.  Covariates are referred to by their names
   in the data frame \code{covar}.} 
 \item{method}{Indicates whether to use multiple imputation or
   Haley-Knott regression.}
 \item{dropone}{If TRUE, do drop-one-term analysis.}
 \item{get.ests}{If TRUE, return estimated QTL effects and their
   estimated variance-covariance matrix.}
 \item{run.checks}{If TRUE, check the input formula and check for
   individuals with missing phenotypes or covariates.}
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

  \bold{The part to get estimated QTL effects is fully working only for
    the case of autosomes in a backcross or intercross.
    In other cases the values returned are based on a
    design matrix that is convenient for calculations but not easily
    interpreted.}

}

\section{Value}{
  An object of class \code{fitqtl}. It may contains as many as three fields:
  \enumerate{
    \item result.full is the ANOVA table as a matrix for the full model
    result. It contains the degree of freedom (df), Sum of squares (SS),
    mean square (MS), LOD score (LOD), percentage of variance explained
    (\%var) and P value (Pvalue).


    \item result.drop is a drop-one-term ANOVA table as a matrix. It
    contains degrees of freedom (df), Type III sum of squares (Type III
    SS), LOD score(LOD), percentage of variance explained (\%var), F
    statistics (F value),  and P values for chi square (Pvalue(chi2))
    and F distribution (Pvalue(F)).

    Note that the degree of freedom, Type III sum of squares, the LOD
    score and the percentage of variance explained are the values
    comparing the full to the sub-model with the term dropped. Also note
    that for imputation method, the percentage of variance explained, the
    the F values and the P values are approximations calculated from the
    LOD score.

    \item ests contains the estimated QTL effects and standard errors.  

  \bold{The part to get estimated QTL effects is fully working only for
    the case of autosomes in a backcross or intercross.
    In other cases the values returned are based on a
    design matrix that is convenient for calculations but not easily
    interpreted.}
    
  }
}

\references{
  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Sen, \'S. and Churchill, G. A. (2001) A statistical framework for quantitative
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
fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp)

# fit model with 3 interacting QTLs interacting
# (performing a drop-one-term analysis)
lod <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1*Q2*Q3)
summary(lod)

\dontrun{
# fit an additive QTL model
lod.add <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1+Q2+Q3)
summary(lod.add)

# fit the model including sex as an interacting covariate
Sex <- data.frame(Sex=fake.f2$pheno$sex)
lod.sex <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1*Q2*Q3*Sex, cov=Sex)
summary(lod.sex)

# fit the same with an additive model
lod.sex.add <- fitqtl(fake.f2, pheno.col=1, qtl, formula=y~Q1+Q2+Q3+Sex, cov=Sex)
summary(lod.sex.add)
}}

\seealso{ \code{\link[qtl]{summary.fitqtl}}, \code{\link[qtl]{makeqtl}},
\code{\link[qtl]{scanqtl}}, \code{\link[qtl]{refineqtl}},
  \code{\link[qtl]{addtoqtl}} ,
  \code{\link[qtl]{dropfromqtl}},
  \code{\link[qtl]{replaceqtl}},
  \code{\link[qtl]{reorderqtl}}}  

\keyword{models}
