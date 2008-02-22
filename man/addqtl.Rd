\name{addqtl}
\alias{addqtl}

\title{Scan for an additional QTL in a multiple-QTL model}

\description{
  Scan for an additional QTL in the context of a
  multiple QTL model.
}

\usage{
addqtl(cross, chr, pheno.col=1, qtl, covar=NULL, formula,
       method=c("imp","hk"), incl.markers=TRUE, verbose=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to be scanned.
   If missing, all chromosomes are scanned.}
 \item{pheno.col}{Column number in the phenotype matrix to be
   used as the phenotype.  One may also give a character string matching
   a phenotype name.}
 \item{qtl}{An object of class \code{qtl}, as output from
   \code{\link[qtl]{makeqtl}}.} 
 \item{covar}{A data.frame of covariates.  These must be strictly
   numeric.} 
 \item{formula}{An object of class \code{\link[stats]{formula}}
   indicating the model to be fitted.  QTLs are referred to as
   \code{Q1}, \code{Q2}, etc.  Covariates are referred to by their names
   in the data frame \code{covar}.  If the new QTL is not included in
   the formula, its main effect is added.}
 \item{method}{Indicates whether to use multiple imputation or
   Haley-Knott regression.}
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.  If \code{\link[qtl]{calc.genoprob}} or
   \code{\link[qtl]{sim.geno}} were run with
   \code{stepwidth="variable"}, we force \code{incl.markers=TRUE}.}
 \item{verbose}{If TRUE, display information about the progress of
   calculations.  If \code{verbose} is an integer > 1, further messages
   from \code{\link[qtl]{scanqtl}} are also displayed.}
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
  \code{\link[qtl]{scanone}} function.  LOD scores are relative to the
  base model (with any terms that include the new QTL omitted).
}

\references{
  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Sen, \'S. and Churchill, G. A. (2001) A statistical framework for quantitative
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

fake.f2 <- sim.geno(fake.f2, n.draws=8, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp)

# scan for an additional QTL
out1 <- addqtl(fake.f2, qtl=qtl, formula=y~Q1+Q2+Q3)
max(out1)

# scan for an additional QTL that interacts with the locus on chr 1
out2 <- addqtl(fake.f2, qtl=qtl, formula=y~Q1*Q4+Q2+Q3)
max(out2)

# plot interaction LOD scores
plot(out2-out1)
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{fitqtl}},
  \code{\link[qtl]{scanqtl}}, \code{\link[qtl]{refineqtl}},
  \code{\link[qtl]{makeqtl}}, \code{\link[qtl]{addtoqtl}},
  \code{\link[qtl]{addpair}}, \code{\link[qtl]{addint}} }

\keyword{models}
