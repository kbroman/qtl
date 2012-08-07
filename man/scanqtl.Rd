\name{scanqtl}
\alias{scanqtl}

\title{General QTL scan}

\description{
  Performs a multiple QTL scan for specified chromosomes and positions
  or intervals, with the possible inclusion of QTL-QTL interactions and/or
  covariates.
}

\usage{
scanqtl(cross, pheno.col=1, chr, pos, covar=NULL, formula,
        method=c("imp","hk"), model=c("normal", "binary"),
        incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
        forceXcovar=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{Column number in the phenotype matrix to be
   used as the phenotype.  One may also give a character string matching
   a phenotype name.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{chr}{Vector indicating the chromosome for each QTL. (These should
   be character strings referring to the chromosomes by name.)}
 \item{pos}{List indicating the positions or intervals on the chromosome
   to be scanned.  Each element should be either a single number (for a
   specific position) or a pair of numbers (for an interval).}
 \item{covar}{A matrix or data.frame of covariates.  These must be
   strictly numeric.}
 \item{formula}{An object of class \code{\link[stats]{formula}}
   indicating the model to be fitted.  (It can also be the character
   string representation of a formula.)  QTLs are indicated as \code{Q1},
   \code{Q2}, etc.  Covariates are indicated by their names in \code{covar}.}
 \item{method}{Indicates whether to use multiple imputation or
   Haley-Knott regression.}
 \item{model}{The phenotype model: the usual model or a model for binary
  traits}
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.  If \code{\link{calc.genoprob}} or
   \code{\link{sim.geno}} were run with
   \code{stepwidth="variable"} or \code{stepwidth="max"}, we force \code{incl.markers=TRUE}.}
 \item{verbose}{If TRUE, give feedback about progress.}
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
  covariate is in involved in an interaction, its main effect are also
  be included.

  Only the interaction terms need to be specifed in the formula. The
  main effects of all input QTLs (as specified by chr and pos) and
  covariates (as specifed by covar) will be included by default. For
  example, if the formula is \code{y~Q1*Q2*Sex}, and there are three
  elements in input \code{chr} and \code{pos} and Sex is one of the
  column names for
  input covariates, the formula used in genome scan will be
  \code{y ~ Q1 + Q2 + Q3 + Sex + Q1:Q2 + Q1:Sex + Q2:Sex + Q1:Q2:Sex}.

  The input \code{pos} is a list or vector to specify the position/range
  of the input chromosomes to be scanned. If it is a vector, it gives the
  precise positions of the QTL on the chromosomes. If it is a list, it will
  contain either the precise positions or a range on the chromosomes. For
  example, consider the case that the input \code{chr = c(1, 6,
    13)}. If \code{pos = c(9.8, 34.0, 18.6)},
  it means to fit a model with QTL on chromosome 1 at 9.8cM, chromosome
  6 at 34cM and  chromosome 13 at 18.6cM.
  If \code{pos = list(c(5,15), c(30,36), 18)}, it
  means to scan chromosome 1 from 5cM to 15cM, chromosome 6 from 30cM to
  36cM, fix the QTL on chromosome 13 at 18cM.
}

\value{
  An object of class \code{scanqtl}. It is a multi-dimensional
  array of LOD scores, with the number of dimension equal to the number
  of QTLs specifed.}

\references{
  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Hao Wu}

\examples{
data(fake.f2)

# take out several QTLs
qc <- c(1, 8, 13)
fake.f2 <- subset(fake.f2, chr=qc)

# imputate genotypes
\dontshow{fake.f2 <- subset(fake.f2, ind=1:50)}
fake.f2 <- calc.genoprob(fake.f2, step=5, err=0.001)

# 2-dimensional genome scan with additive 3-QTL model
pos <- list(c(15,35), c(45,65), 28)
result <- scanqtl(fake.f2, pheno.col=1, chr=qc, pos=pos,
                  formula=y~Q1+Q2+Q3, method="hk")

# image of the results
# chr locations
chr1 <- as.numeric(matrix(unlist(strsplit(colnames(result),"@")),
                   ncol=2,byrow=TRUE)[,2])
chr8 <- as.numeric(matrix(unlist(strsplit(rownames(result),"@")),
                   ncol=2,byrow=TRUE)[,2])
# image plot
image(chr1, chr8, t(result), las=1, col=rev(rainbow(256,start=0,end=2/3)))

# do the same, allowing the QTLs on chr 1 and 13 to interact
result2 <- scanqtl(fake.f2, pheno.col=1, chr=qc, pos=pos,
                   formula=y~Q1+Q2+Q3+Q1:Q3, method="hk")
# image plot
image(chr1, chr8, t(result2), las=1, col=rev(rainbow(256,start=0,end=2/3)))
}

\seealso{ \code{\link{fitqtl}}, \code{\link{makeqtl}},
  \code{\link{refineqtl}} }

\keyword{models}
