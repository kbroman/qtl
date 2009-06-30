\name{scantwo}
\alias{scantwo}

\title{Two-dimensional genome scan with a two-QTL model}

\description{
  Perform a two-dimensional genome scan with a two-QTL model, with
  possible allowance for covariates. 
}

\usage{
scantwo(cross, chr, pheno.col=1, model=c("normal","binary"),
        method=c("em","imp","hk","mr","mr-imp","mr-argmax"),
        addcovar=NULL, intcovar=NULL, weights=NULL,
        use=c("all.obs", "complete.obs"), 
        incl.markers=FALSE, clean.output=FALSE, 
        clean.nmar=1, clean.distance=0,
        maxit=4000, tol=1e-4,
        verbose=TRUE, n.perm, perm.strata=NULL,
        assumeCondIndep=FALSE, batchsize=250, n.cluster=1)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes for which LOD
   scores should be calculated.  This should be a vector of character
   strings referring to chromosomes by name; numeric values are
   converted to strings.  Refer to chromosomes with a preceding \code{-}
   to have all chromosomes but those considered.  A logical (TRUE/FALSE)
   vector may also be used.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This can be a vector of integers; for methods
   \code{"hk"} and \code{"imp"} this can be considerably faster than doing
   them one at a time.  One may also give character strings matching
   the phenotype names.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{model}{The phenotypic model: the usual normal model or a model for
   binary traits.}
 \item{method}{Indicates whether to use the
   the EM algorithm, imputation, Haley-Knott regression, or marker
   regression.  Marker regression is performed either by dropping
   individuals with missing genotypes (\code{"mr"}), or by first filling
   in missing data using a single imputation (\code{"mr-imp"}) or by the
   Viterbi algorithm (\code{"mr-argmax"}).}
 \item{addcovar}{Additive covariates.}
 \item{intcovar}{Interactive covariates (interact with QTL genotype).}
 \item{weights}{Optional weights of individuals.  Should be either NULL
   or a vector of length n.ind containing positive weights. Used only
   in the case \code{model="normal"}.}
 \item{use}{In the case that multiple phenotypes are selected to be
   scanned, this argument indicates whether to use all individuals, 
   including those missing some phenotypes, or just those individuals
   that have data on all selected phenotypes.}
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.  If \code{\link[qtl]{calc.genoprob}} or
   \code{\link[qtl]{sim.geno}} were run with
   \code{stepwidth="variable"}, we force \code{incl.markers=TRUE}.}
 \item{clean.output}{If TRUE, clean the output with
   \code{\link[qtl]{clean.scantwo}}, replacing LOD scores for pairs of
   positions that are not well separated with 0.  In permutations, this
   will be done for each permutation replicate.  This can be important
   for the case of \code{method="em"}, as there can be difficulty with
   algorithm convergence in these regions.}
 \item{clean.nmar}{If \code{clean.output=TRUE}, this is the number of
 markers that must separate two positions.}
 \item{clean.distance}{If \code{clean.output=TRUE}, this is the cM distance
 that must separate two positions.}
 \item{maxit}{Maximum number of iterations; used
   only with method \code{"em"}.}
 \item{tol}{Tolerance value for determining convergence; used only with
   method \code{"em"}.} 
 \item{verbose}{If TRUE, display information about the progress of
   calculations.  For method \code{"em"}, if \code{verbose} is an integer
   above 1, further details on the progress of the algorithm will be
   displayed.}
 \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
 \item{perm.strata}{If \code{n.perm} > 0, this may be used to perform a
   stratified permutation test.  This should be a vector with the same
   number of individuals as in the cross data.  Unique values indicate
   the individual strata, and permutations will be performed within the
   strata.}
 \item{assumeCondIndep}{If TRUE, assume conditional independence of QTL
   genotypes given marker genotypes.  This is an approximation, but it
   may speed things up.}
 \item{batchsize}{The number of phenotypes (or permutations) to be run
   as a batch; used only for methods \code{"hk"} and \code{"imp"}.}
 \item{n.cluster}{If the package \code{\link[snow:00Index]{snow}} is available and
   \code{n.perm} > 0, permutations are run in parallel using this number
   of nodes.}
}

\details{
  Standard interval mapping (\code{method="em"}) and Haley-Knott
  regression (\code{method="hk"}) require that multipoint genotype probabilities are
  first calculated using \code{\link[qtl]{calc.genoprob}}.  The
  imputation method uses the results of \code{\link[qtl]{sim.geno}}. 

  The method \code{"em"} is standard interval mapping by the EM algorithm
  (Dempster et al. 1977; Lander and Botstein 1989).  Marker regression
  (\code{method="mr"}) is simply linear regression of phenotypes on
  marker genotypes  (individuals with missing genotypes are
  discarded). Haley-Knott regression (\code{method="hk"}) uses the
  regression of phenotypes on multipoint genotype probabilities.  The
  imputation method (\code{method="imp"}) uses the pseudomarker 
  algorithm described by Sen and Churchill (2001).
  
  Individuals with missing phenotypes are dropped.

  In the presence of covariates, the full model is 
  \deqn{y = \mu + \beta_{q_1} + \beta_{q_2} + \beta_{q_1 \times q_2} +
    A \gamma + Z \delta_{q_1} + Z \delta_{q_2} + Z \delta_{q_1 \times q_2} +
    \epsilon}{y = m + b[q1] + b[q2] + b[q1 x q2] + A g + Z d[q1] + Z
    d[q2] + Z d[q1 x q2] + e}
  where \eqn{q_1}{q1} and \eqn{q_2}{q2} are the unknown QTL genotypes at two
  locations, \emph{A} is a matrix of covariates, and \emph{Z} is a
  matrix of covariates that interact with QTL genotypes.  The columns of
  \emph{Z} are forced to be contained in the matrix \emph{A}.

  The above full model is compared to the additive QTL model, 
  \deqn{y = \mu + \beta_{q_1} + \beta_{q_2} + A \gamma + Z \delta_{q_1} + Z
    \delta_{q_2} + \epsilon}{y = m + b[q1] + b[q2] + A g + Z d[q1] + Z
    d[q2] + e}
  and also to the null model, with no QTL, 
  \deqn{y = \mu + A \gamma + \epsilon}{y = m + A g + e}

  In the case that \code{n.perm} is specified, the R function
  \code{scantwo} is called repeatedly.
}

\section{X chromosome}{
The X chromosome must be treated specially in QTL mapping.

As in \code{\link[qtl]{scanone}}, if both males and females are
included, male hemizygotes are allowed to be different from female
homozygotes, and the null hypothesis must be changed in order to ensure
that sex- or pgm-differences in the phenotype do not results in spurious
linkage to the X chromosome.  (See the help file for
\code{\link[qtl]{scanone}}.) 
}

\value{
  If \code{n.perm} is missing, the function returns a list with class
  \code{"scantwo"} and containing three components.  The first component
  is a matrix of dimension [tot.pos x tot.pos]; the upper triangle
  contains the LOD scores for the additive model, and the lower triangle
  contains the LOD scores for the full model.  The diagonal contains the
  results of \code{\link[qtl]{scanone}}. The second component of the
  output is a data.frame indicating the locations at which the two-QTL
  LOD scores were calculated.  The first column is the chromosome
  identifier, the second column is the position in cM, the third column
  is a 1/0 indicator for ease in later pulling out only the equally
  spaced positions, and the fourth column indicates whether the position
  is on the X chromosome or not.  The final component is a version of
  the results of \code{\link[qtl]{scanone}} including sex and/or cross
  direction as additive covariates, which is needed for a proper
  calculation of conditional LOD scores. 

  If \code{n.perm} is specified, the function returns a list with six
  different LOD scores from each of the permutation replicates. 
  First, the maximum LOD score for the full model (two QTLs plus an
  interaction).  Second, for each pair of
  chromosomes, we take the difference between the full LOD and the
  maximum single-QTL LOD for those two chromosomes, and then maximize
  this across chromosome pairs.  Third, for each pair of chromosomes we
  take the difference between the maximum full LOD and the maximum
  additive LOD, and then maximize this across chromosome pairs.  Fourth,
  the maximum LOD score for the additive QTL model.  Fifth, for each
  pair of chromosomes, we take the difference between the additive LOD
  and the maximum single-QTL LOD for those two chromosomes, and then
  maximize this across chromosome pairs.  Finally, the maximum
  single-QTL LOD score (that is, from a single-QTL scan).  The latter is
  not used in \code{\link[qtl]{summary.scantwo}}, but does get
  calculated at each permutation, so we include it for the sake of
  completeness. 
}

\references{
  Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for
  quantitative trait mapping.  \emph{Genetics} \bold{138}, 963--971.

  Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum
  likelihood from incomplete data via the EM algorithm.  \emph{J. Roy.
    Statist. Soc.} B, \bold{39}, 1--38.

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Lander, E. S. and Botstein, D. (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps.  \emph{Genetics}
  \bold{121}, 185--199.

  Sen, \'S. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.

  Soller, M., Brody, T. and Genizi, A. (1976) On the power of experimental
  designs for the detection of linkage between marker loci and
  quantitative loci in crosses between inbred lines.
  \emph{Theor. Appl. Genet.} \bold{47}, 35--39. 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}; Hao Wu}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
fake.f2 <- calc.genoprob(fake.f2, step=5)
out.2dim <- scantwo(fake.f2, method="hk")
plot(out.2dim)

# permutations
\dontshow{permo.2dim <- scantwo(fake.f2, method="hk", n.perm=2)}
\dontrun{permo.2dim <- scantwo(fake.f2, method="hk", n.perm=1000)}
summary(permo.2dim, alpha=0.05)

# summary with p-values
summary(out.2dim, perms=permo.2dim, pvalues=TRUE,
        alphas=c(0.05, 0.10, 0.10, 0.05, 0.10))

# covariates
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=16:17)}
fake.bc <- calc.genoprob(fake.bc, step=10)

ac <- pull.pheno(fake.bc, c("sex","age"))
ic <- pull.pheno(fake.bc, "sex")

out <- scantwo(fake.bc, method="hk", pheno.col=1,
               addcovar=ac, intcovar=ic)
plot(out)
}

\seealso{ \code{\link[qtl]{plot.scantwo}}, \code{\link[qtl]{summary.scantwo}},
  \code{\link[qtl]{scanone}}, \code{\link[qtl]{max.scantwo}},
  \code{\link[qtl]{summary.scantwoperm}},
  \code{\link[qtl]{c.scantwoperm}} }

\keyword{models}
