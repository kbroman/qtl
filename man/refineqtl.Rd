\name{refineqtl}
\alias{refineqtl}

\title{Refine the positions of QTL}

\description{
  Iteratively scan the positions for QTL in the context of a multiple
  QTL model, to try to identify the positions with maximum likelihood,
  for a fixed QTL model.
}

\usage{
refineqtl(cross, pheno.col=1, qtl, chr, pos, qtl.name, covar=NULL, formula,
          method=c("imp","hk"), model=c("normal", "binary"), verbose=TRUE, maxit=10,
          incl.markers=TRUE, keeplodprofile=TRUE, tol=1e-4,
          maxit.fitqtl=1000, forceXcovar=FALSE)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{Column number in the phenotype matrix to be
   used as the phenotype.  One may also give a character string matching
   the phenotype name.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{qtl}{A QTL object, as produced by \code{\link{makeqtl}},
   containing the positions of the QTL.  Provide either \code{qtl} or
   the pair \code{chr} and \code{pos}.}
 \item{chr}{Vector indicating the chromosome for each QTL; if \code{qtl}
 is provided, this should not be.}
 \item{pos}{Vector indicating the positions for each QTL; if \code{qtl}
   is provided, this should not be.}
 \item{qtl.name}{Optional user-specified name for each QTL.  If
   \code{qtl} is provided, this should not be.}
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
 \item{verbose}{If TRUE, give feedback about progress.  If
   \code{verbose} is an integer > 1, further messages from
   \code{\link{scanqtl}} are also displayed.}
 \item{maxit}{Maximum number of iterations.}
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.}
 \item{keeplodprofile}{If TRUE, keep the LOD profiles from the last
   iteration as attributes to the output.}
 \item{tol}{Tolerance for convergence for the binary trait model.}
 \item{maxit.fitqtl}{Maximum number of iterations for fitting the binary trait
   model.}
 \item{forceXcovar}{If TRUE, force inclusion of X-chr-related covariates
   (like sex and cross direction).}
}

\details{
  QTL positions are optimized, within the context of a fixed QTL model,
  by a scheme described in Zeng et al. (1999).  Each QTL is considered
  one at a time (in a random order), and a scan is performed, allowing
  the QTL to vary across its chromosome, keeping the positions of all
  other QTL fixed.  If there is another QTL on the chromosome, the
  position of the floating QTL is scanned from the end of the chromosome
  to the position of the flanking QTL.  If the floating QTL is between
  two QTL on a chromosome, its position is scanned between those two QTL
  positions.  Each QTL is moved to the position giving the highest
  likelihood, and the entire process is repeated until no further
  improvement in likelihood can be obtained.

  One may provide either a \code{qtl} object (as produced by
  \code{\link{makeqtl}}), or vectors \code{chr} and \code{pos}
  (and, optionally, \code{qtl.name}) indicating the positions of the
  QTL.

  If a \code{qtl} object is provided, QTL that do not appear in
  the model \code{formula} are ignored, but they remain part of the QTL
  object that is output.
}

\value{
  An object of class \code{qtl}, with QTL placed in their new positions.

  If \code{keeplodprofile=TRUE}, LOD profiles from the last pass through
  the refinement algorithm are retained as an attribute,
  \code{"lodprofile"}, to the object.  These may be plotted with
  \code{\link{plotLodProfile}}.
}

\references{
  Zeng, Z.-B., Kao, C.-H., and Basten, C. J. (1999) Estimating the
  genetic architecture of quantitative traits.  \emph{Genet. Res.}
  \bold{74}, 279--289.

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=c(2,5))}
fake.bc <- calc.genoprob(fake.bc, step=2)
qtl <- makeqtl(fake.bc, chr=c(2,5), pos=c(32.5, 17.5), what="prob")
rqtl <- refineqtl(fake.bc, qtl=qtl, method="hk")
}

\seealso{ \code{\link{fitqtl}}, \code{\link{makeqtl}},
  \code{\link{scanqtl}}, \code{\link{addtoqtl}},
  \code{\link{dropfromqtl}}, \code{\link{replaceqtl}},
  \code{\link{plotLodProfile}} }

\keyword{utilities}
