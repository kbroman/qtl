\name{addpair}
\alias{addpair}

\title{Scan for an additional pair of QTL in a multiple-QTL model}

\description{
  Scan for an additional pair of QTL in the context of a
  multiple QTL model.
}

\usage{
addpair(cross, chr, pheno.col=1, qtl, covar=NULL, formula,
        method=c("imp","hk"), model=c("normal", "binary"),
        incl.markers=FALSE, verbose=TRUE, tol=1e-4, maxit=1000,
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
   in the data frame \code{covar}.  If the new QTL are not included in
   the formula, a two-dimensional scan as in \code{\link{scantwo}}
   is performed.}
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

  If neither of the two new QTL are indicated in the \code{formula}, we
  perform a two-dimensional scan as in \code{\link{scantwo}}.  That
  is, for each pair of QTL positions, we fit two models: two
  additive QTL added to the \code{formula}, and two interacting QTL
  added to the \code{formula}.

  If the both of the new QTL are indicated in the \code{formula}, that
  particular model is fit, with the positions of the new QTL allowed to
  vary across the genome.  If just one of the QTL is indicated in the
  \code{formula}, a main effect for the other is added, and that
  particular model is fit, again with the positions of both QTL varying.
  Note that in this case the LOD scores are not analogous to those
  produced by \code{\link{scantwo}}.  Thus, there slightly modified forms
  for the plots (produced by \code{\link{plot.scantwo}}) and
  summaries (produced by \code{\link{summary.scantwo}} and
  \code{\link{max.scantwo}}).  In the plot, the x-axis is to be
  interpreted as the position of the first of the new QTL, and the
  y-axis is to be interpreted as the position of the second of the new
  QTL.  In the summaries, we give the single best pair of positions on
  each pair of chromosomes, and give LOD scores comparing that pair of
  positions to the base model (without each of these QTL), and to the
  base model plus one additional QTL on one or the other of the chromosomes.
}

\section{Value}{
  An object of class \code{scantwo}, as produced by
  \code{\link{scantwo}}.

  If neither of the new QTL were indicated
  in the \code{formula}, the result is just as in
  \code{\link{scantwo}}, though with LOD scores relative to the
  base model (omitting the new QTL).

  Otherwise, the results are contained in what would ordinarily be in
  the full and additive LOD scores, with the additive LOD scores
  corresponding to the case that the first of the new QTL is to the left
  of the second of the new QTL, and the full LOD scores corresponding to
  the case that the first of the new QTL is to the right of the second
  of the new QTL.  Because the structure of the LOD scores in this
  case is different from those output by \code{\link{scantwo}}, we
  include, in this case, an attribute \code{"addpair"=TRUE}.  (We also
  require results of single-dimensional scans, omitting each of the two
  new QTL from the formula, one at a time; these are included as
  attributes \code{"lod.minus1"} and \code{"lod.minus2"}.)  The
  results are then treated somewhat differently by
  \code{\link{summary.scantwo}}, \code{\link{max.scantwo}},
  and \code{\link{plot.scantwo}}.  See the Details section.
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
# A totally contrived example to show some of what you can do

# simulate backcross data with 3 chromosomes (names "17", "18", "19")
#   one QTL on chr 17 at 40 cM
#   one QTL on chr 18 at 30 cM
#   two QTL on chr 19, at 10 and 40 cM
data(map10)
model <- rbind(c(1,40,0), c(2,30,0), c(3,10,0), c(3,40,0))
\dontrun{fakebc <- sim.cross(map10[17:19], model=model, type="bc", n.ind=250)}
\dontshow{fakebc <- sim.cross(map10[17:19], model=model, type="bc", n.ind=25)}

# het at QTL on 17 and 1st QTL on 19 increases phenotype by 1 unit
# het at QTL on 18 and 2nd QTL on 19 decreases phenotype by 1 unit
qtlgeno <- fakebc$qtlgeno
phe <- rnorm(nind(fakebc))
w <- qtlgeno[,1]==2 & qtlgeno[,3]==2
phe[w] <- phe[w] + 1
w <- qtlgeno[,2]==2 & qtlgeno[,4]==2
phe[w] <- phe[w] - 1
fakebc$pheno[,1] <- phe

\dontrun{fakebc <- calc.genoprob(fakebc, step=2, err=0.001)
}\dontshow{fakebc <- calc.genoprob(fakebc, step=0, err=0.001)}

# base model has QTLs on chr 17 and 18
qtl <- makeqtl(fakebc, chr=c("17", "18"), pos=c(40,30), what="prob")

# scan for an additional pair of QTL, one interacting with the locus
#     on 17 and one interacting with the locus on 18
out.ap <- addpair(fakebc, qtl=qtl, formula = y~Q1*Q3 + Q2*Q4, method="hk")

max(out.ap)
summary(out.ap)
plot(out.ap)
}

\seealso{ \code{\link{addint}}, \code{\link{addqtl}},
  \code{\link{fitqtl}}, \code{\link{makeqtl}},
  \code{\link{scanqtl}}, \code{\link{refineqtl}},
  \code{\link{makeqtl}}, \code{\link{scantwo}},
  \code{\link{addtoqtl}} }

\keyword{models}
