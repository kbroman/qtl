\name{stepwiseqtl}
\alias{stepwiseqtl}

\title{Stepwise selection for multiple QTL}

\description{
  Performs forward/backward selection to identify a multiple QTL model,
  with model choice made via a penalized LOD score, with separate
  penalties on main effects and interactions.
}

\usage{
stepwiseqtl(cross, chr, pheno.col=1, qtl, formula, max.qtl=10, covar=NULL,
            method=c("imp", "hk"), model=c("normal", "binary"),
            incl.markers=TRUE, refine.locations=TRUE,
            additive.only=FALSE, scan.pairs=FALSE, penalties,
            keeplodprofile=FALSE, keeptrace=FALSE, verbose=TRUE,
            tol=1e-4, maxit=1000)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider in
   search for QTL.  This should be a vector of character strings
   referring to chromosomes by name; numeric values are converted to
   strings.  Refer to chromosomes with a preceding \code{-} to have all
   chromosomes but those considered.  A logical (TRUE/FALSE) vector may
   also be used.} 
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  One may also give character strings matching
   the phenotype names.  Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{qtl}{Optional QTL object (of class \code{"qtl"}, as created by
   \code{\link{makeqtl}}) to use as a starting point.}
 \item{formula}{Optional formula to define the QTL model to be used as a
   starting point.} 
 \item{max.qtl}{Maximum number of QTL to which forward selection should
   proceed.}
 \item{covar}{Data frame of additive covariates.}
 \item{method}{Indicates whether to use multiple imputation or
   Haley-Knott regression.}
 \item{model}{The phenotype model: the usual model or a model for binary
  traits} 
 \item{incl.markers}{If FALSE, do calculations only at points on an
   evenly spaced grid.}
 \item{refine.locations}{If TRUE, use \code{\link{refineqtl}} to
   refine the QTL locations after each step of forward and backward
   selection.}
 \item{additive.only}{If TRUE, allow only additive QTL models; if FALSE,
 consider also pairwise interactions among QTL.}
 \item{scan.pairs}{If TRUE, perform a two-dimensional, two-QTL scan at
   each step of forward selection.}
 \item{penalties}{Vector of three values indicating the penalty on main
   effects and heavy and light penalties on interactions.  See the
   Details below. If missing, default values are used that are based on
   simulations of backcrosses and intercrosses with genomes modeled
   after that of the mouse.}
 \item{keeplodprofile}{If TRUE, keep the LOD profiles from the last
   iteration as attributes to the output.}
 \item{keeptrace}{If TRUE, keep information on the sequence of models
   visited through the course of forward and backward selection as an
   attribute to the output.}
 \item{verbose}{If TRUE, give feedback about progress.  If
   \code{verbose} is an integer > 1, even more information is printed.} 
 \item{tol}{Tolerance for convergence for the binary trait model.}
 \item{maxit}{Maximum number of iterations for fitting the binary trait
   model.} 
}

\section{Methods}{
  \bold{\code{imp}}: multiple imputation is used, as described by Sen
  and Churchill (2001).

  \bold{\code{hk}}: Haley-Knott regression is used (regression of the
  phenotypes on the multipoint QTL genotype probabilities), as described
  by Haley and Knott (1992).
}

\details{
  We seek to identify the model with maximal penalized LOD score.  The
  penalized LOD score, defined in Manichaikul et al. (2009),
  is the LOD score for the model (the \eqn{\log_{10}}{log10} likelihood
  ratio comparing the model to the null model with no QTL) with
  penalties on the number of QTL and QTL:QTL interactions.

  We consider QTL models allowing pairwise interactions among QTL but
  with an enforced hierarchy in which inclusion of a pairwise
  interaction requires the inclusion of both of the corresponding main
  effects.  Additive covariates may be included, but currently we do not
  explore QTL:covariate interactions.  Also, the penalized LOD score
  criterion is currently defined only for autosomal loci, and results
  with the X chromosome should be considered with caution.

  The penalized LOD score is of the form \eqn{pLOD(\gamma) = LOD(\gamma)
  - T_m p_m - T_h p_h - T_l p_l}{pLOD(g) = LOD(g)
  - Tm pm - Th ph - Tl pl} where \eqn{\gamma}{g} denotes a model,
  \eqn{p_m}{pm} is the number of QTL in the model ("main effects"),
  \eqn{p_h}{ph} is the number of pairwise interactions that will be
  given a heavy interaction penalty, \eqn{p_l}{pl} is the number of pairwise
  interactions that will be given a light interaction penalty,
  \eqn{T_m}{Tm} is the penalty on main effects, \eqn{T_h}{Th} is the heavy
  interaction penalty, and \eqn{T_l}{Tl} is the light interaction
  penalty.  The \code{penalties} argument is the vector \eqn{(T_m, T_h,
  T_l)}{(Tm, Th, Tl)}.  If \eqn{T_l}{Tl} is missing (\code{penalties}
  has a vector of length 2), we assume \eqn{T_l = T_h}{Tl = Th}, and so
  all pairwise interactions are assigned the same penalty.  

  The "heavy" and "light" interaction penalties can be a bit
  confusing.  Consider the clusters of QTL that are connected via one or
  more pairwise interactions.  To each such cluster, we assign at most
  one "light" interaction penalty, and give all other pairwise
  interactions the heavy interaction penalty.  In other words, if 
  \eqn{p_i}{pi} is the total number of pairwise interactions for a QTL
  model, we let \eqn{p_l}{pl} be the number of clusters of connected QTL
  with at least one pairwise interaction, and then let \eqn{p_h - p_i -
    p_l}{ph = pi - pl}.

  Let us give an explicit example.  Consider a model with 6 QTL, and
  with interactions between QTL 2 and 3, QTL 4 and 5 and QTL 4 and 6
  (so we have the model formula
  \code{y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q2:Q3 + Q4:Q5 + Q4:Q6}).
  There are three clusters of connected QTL: (1), (2,3) and (4,5,6).  We
  would assign 6 main effect penalties (\eqn{T_m}{Tm}), 2 light
  interaction penalties (\eqn{T_l}{Tl}), and 1 heavy interaction penalty
  (\eqn{T_h}{Th}).  
  
  Manichaikul et al. (2009) described a system for deriving the
  three penalties on the basis of permutation results from a
  two-dimensional, two-QTL genome scan (as calculated with
  \code{\link{scantwo}}).  These may be calculated with the
  function \code{\link{calc.penalties}}.

  A forward/backward search method is used, with the aim to optimize the
  penalized LOD score criterion.  That is, we seek to identify the
  model with maximal the penalized LOD score.  The search algorithm was
  based closely on an algorithm described by Zeng et al. (1999).

  We use forward selection to a model of moderate size (say 10 QTL),
  followed by backward elimination all the way to the null model.  The
  chosen model is that which optimizes the penalized LOD score
  criterion, among all models visited.  The detailed algorithm is as
  follows.  Note that if \code{additive.only=TRUE}, no pairwise
  interactions are considered.

  \enumerate{
    \item Start at the null model, and perform a single-QTL genome scan,
    and choose the position giving the largest LOD score.  If
    \code{scan.pairs=TRUE}, start with a two-dimensional, two-QTL genome
    scan instead.  If an initial QTL model were defined through the
    arguments \code{qtl} and \code{formula}, start with this model and
    jump immediately to step 2.

    \item With a fixed QTL model in hand:
    \enumerate{
      \item Scan for an additional additive QTL.

      \item For each QTL in the current model, scan for an additional
      interacting QTL.

      \item If there are \eqn{\ge}{>=} 2 QTL in the current model,
      consider adding one of the possible pairwise interactions.

      \item If \code{scan.pairs=TRUE} perform a two-dimensional, two-QTL
      scan, seeking to add a pair of novel QTL, either additive or
      interacting.

      \item Step to the model that gives the largest value for the model
      comparison criterion, among those considered at the current step.
      }

      \item Refine the locations of the QTL in the current model (if
      \code{refine.locations=TRUE}).  

      \item Repeat steps 2 and 3 up to a model with some pre-determined
      number of loci.

      \item Perform backward elimination, all the way back to the null
      model.  At each step, consider dropping one of the current main
      effects or interactions; move to the model that maximizes the model
      comparison criterion, among those considered at this step.  Follow
      this with a refinement of the locations of the QTL.

      \item Finally, choose the model having the largest model comparison
      criterion, among all models visited.
  }

  In this forward/backward algorithm, it is likely best to build up to
  an overly large model and then prune it back.  Note that there is no
  "stopping rule"; the chosen model is that which optimizes the model
  comparison criterion, among all models visited.  The search can be
  time consuming, particularly if a two-dimensional scan is performed at
  each forward step.  Such two-dimensional scans may be useful for
  identifying QTL linked in repulsion (having effects of opposite sign)
  or interacting QTL with limited marginal effects, but our limited
  experience suggests that they are not necessary; important linked or
  interacting QTL pairs can be picked up in the forward selection to a
  large model, and will be retained in the backward elimination phase.
}

\value{
  The output is a representation of the best model, as measured by the
  penalized LOD score (see Details), among all models visited.
  This is QTL object (of class \code{"qtl"}, as produced by
  \code{\link{makeqtl}}), with attributes \code{"formula"},
  indicating the model formula, and \code{"pLOD"} indicating the
  penalized LOD score.

  If \code{keeplodprofile=TRUE}, LOD profiles from the last pass through
  the refinement algorithm are retained as an attribute,
  \code{"lodprofile"}, to the object.  These may be plotted with
  \code{\link{plotLodProfile}}.  

  If \code{keeptrace=TRUE}, the output will contain an attribute
  \code{"trace"} containing information on the best model at each step
  of forward and backward elimination.  This is a list of objects of
  class \code{"compactqtl"}, which is similar to a QTL object (as
  produced by \code{\link{makeqtl}}) but containing just 
  a vector of chromosome IDs and positions for the QTL.  Each will also
  have attributes \code{"formula"} (containing the model formula) and
  \code{"pLOD"} (containing the penalized LOD score.
}

\references{
  Manichaikul, A., Moon, J. Y., Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}, Yandell, B. S. and Broman,
  K. W. (2009) A model selection approach for the identification of
  quantitative trait loci in experimental crosses, allowing
  epistasis. \emph{Genetics}, \bold{181}, 1077--1086.

  Broman,  K. W. and Speed, T. P. (2002) A model selection approach for the
  identification of quantitative trait loci in experimental crosses (with
  discussion). \emph{J Roy Stat Soc B} \bold{64}, 641--656, 731--775.

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.

  Zeng, Z.-B., Kao, C.-H. and Basten, C. J. (1999) Estimating the
  genetic architecture of quantitative traits.  \emph{Genetical
  Research}, \bold{74}, 279--289.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=c(1,2,3,5))}
\dontrun{fake.bc <- calc.genoprob(fake.bc, step=2.5)
}\dontshow{fake.bc <- calc.genoprob(fake.bc, step=0)}
outsw <- stepwiseqtl(fake.bc, max.qtl=3, method="hk", keeptrace=TRUE)

# best model
outsw
plotModel(outsw)

# path through model space
thetrace <- attr(outsw, "trace")

# plot of these
par(mfrow=c(3,3))
for(i in seq(along=thetrace))
  plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))
}

\seealso{ \code{\link{calc.penalties}},
  \code{\link{plotModel}}, \code{\link{makeqtl}},
  \code{\link{fitqtl}}, \code{\link{refineqtl}},
  \code{\link{addqtl}}, \code{\link{addpair}} }

\keyword{models}
