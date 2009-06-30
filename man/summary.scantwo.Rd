\name{summary.scantwo}
\alias{summary.scantwo}

\title{Summarize the results of a two-dimensional genome scan}

\description{
  Summarize the interesting aspects of the results of \code{\link[qtl]{scantwo}}.
}

\usage{
\method{summary}{scantwo}(object, thresholds, 
        what=c("best", "full", "add", "int"),
        perms, alphas, lodcolumn=1, pvalues=FALSE,
        df=FALSE, allpairs=TRUE, \dots)
}

\arguments{
 \item{object}{An object of class \code{scantwo}, the output of the
   function \code{\link[qtl]{scantwo}}.}
 \item{thresholds}{A vector of length 5, giving LOD thresholds
   for the full, conditional-interactive, interaction, additive, and
   conditional-additive LOD scores. See Details, below.}
 \item{what}{Indicates for which LOD score the maximum should be
   reported. See Details, below.} 
 \item{perms}{Optional permutation results used to derive thresholds or
   to calculate genome-scan-adjusted p-values.  This must be consistent
   with the \code{object} input, in that it must have the same number of
   LOD score columns, though it can have just one column of permutation
   results, in which case they are assumed to apply to any chosen LOD
   score column.}
 \item{alphas}{If perms are included, these are the significance levels used
   to calculate thresholds for determining which peaks to pull out.  It
   should be a vector of length 5, giving significance levels 
   for the full, conditional-interactive, interaction, additive, and
   conditional-additive LOD scores.  (It can also be a single number, in
   which case it is assumed that the same value is used for all five LOD
   scores.)  If \code{thresholds} is specified, \code{alphas} should not
   be.} 
 \item{lodcolumn}{If the scantwo results contain LOD scores for
    multiple phenotypes, this argument indicates which to use in the
    summary.  Only one LOD score column may be considered at a time.} 
 \item{pvalues}{If TRUE, include columns with genome-scan-adjusted
   p-values in the results.  This requires that \code{perms} be
   provided.} 
 \item{df}{If TRUE, the degrees of freedom associated with the LOD
   scores are shown.} 
 \item{allpairs}{If TRUE, all pairs of chromosomes are considered.  If
   FALSE, only self-self pairs are considered, so that one may more
   conveniently check for possible linked QTL.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.scantwo}, to be printed by
  \code{print.summary.scantwo};
}

\details{
  If \code{what="best"}, we calculate, for each pair of chromosomes, the
  maximum LOD score for the full model (two QTL plus interaction) and
  the maximum LOD score for the additive model.  The difference between
  these is a LOD score for a test for interaction.  We also calculate
  the difference between the maximum full LOD and the maximum single-QTL LOD
  score for the two chromosomes; this is the LOD score for a test for a
  second QTL, allowing for epistasis, which we call either the
  conditional-interactive or "fv1" LOD score.  Finally,
  we calculate the difference between the maximum additive LOD score and
  the maximum single-QTL LOD score for the two chromosomes; this is the
  LOD score for a test for a second QTL, assuming that the two QTL act
  additively, which we call either the conditional-additive or "av1" LOD
  score.  Note that the maximum full LOD and additive LOD are allowed to
  occur in different places. 

  If \code{what="full"}, we find the maximum full LOD and extract the
  additive LOD at the corresponding pair of positions; we derive
  the other three LOD scores for that fixed pair of positions.

  If \code{what="add"}, we find the maximum additive LOD and extract the
  full LOD  at the corresponding pair of positions; we derive
  the other three LOD scores for that fixed pair of positions.

  If \code{what="int"}, we find the pair of positions for which the
  difference between the full and additive LOD scores is largest, and
  then calculate the five LOD scores at that pair of positions.

  If \code{thresholds} or \code{alphas} is provided (and note that when
  \code{alphas} is provided, \code{perms} must also), we extract just
  those pairs of chromosomes for which either (a) the full LOD score
  exceeds its thresholds and either the conditional-interactive LOD or
  the interaction LOD exceed their threshold, or (b) the additive LOD
  score exceeds its threshold and the conditional-additive LOD exceeds
  its threshold.  The thresholds or alphas must be given in the order
  full, cond-int, int, add, cond-add. 

  Thresholds may be obtained by a permutation test with
  \code{\link[qtl]{scantwo}}, but these are extremely time-consuming.
  For a mouse backcross, we suggest the thresholds (6.0, 4.7, 4.4, 4.7,
  2.6) for the full, conditional-interactive, interaction, additive, and
  conditional-additive LOD scores, respectively.  
  For a mouse intercross, we suggest the thresholds (9.1, 7.1, 6.3, 6.3,
  3.3) for the full, conditional-interactive, interaction, additive, and
  conditional-additive LOD scores, respectively.  These were obtained by
  10,000 simulations of crosses with 250 individuals, markers at a 10 cM
  spacing, and analysis by Haley-Knott regression.
}

\section{Output of addpair}{

  \bold{Note} that, for output from \code{\link[qtl]{addpair}} in which the
  new loci are indicated explicitly in the formula, the summary provided
  by \code{summary.scantwo} is somewhat special.

  All arguments except \code{allpairs} and \code{thresholds} (and, of
  course, the input \code{object}) are ignored.

  If the formula is symmetric in the two new QTL, the output has just two LOD
  score columns: \code{lod.2v0} comparing the full model to the model
  with neither of the new QTL, and \code{lod.2v1} comparing the full
  model to the model with just one new QTL.

  If the formula is \emph{not} symmetric in the two new QTL, the output
  has three LOD score columns: \code{lod.2v0} comparing the full model
  to the model with neither of the new QTL, \code{lod.2v1b} comparing
  the full model to the model in which the first of the new QTL is
  omitted, and \code{lod.2v1a} comparing the full model to the model
  with the second of the new QTL omitted.  

  The \code{thresholds} argument should have length 1 or 2, rather than
  the usual 5.  Rows will be retained if \code{lod.2v0} is greater than
  \code{thresholds[1]} and \code{lod.2v1} (or either of \code{lod.2v1a}
  or \code{lod.2v1b}) is greater than \code{thresholds[2]}.  (If a
  single thresholds is given, we assume that \code{thresholds[2]==0}.)
}

\section{The older version}{
  The previous version of this function is still available, though it is
  now named \code{\link[qtl]{summary.scantwo.old}}.

  We much prefer the revised function.  However, while we are confident
  that this function (and the permutations in
  \code{\link[qtl]{scantwo}}) are calculating the relevant statistics,
  the appropriate significance levels for these relatively complex
  series of statistical tests is not yet completely clear.
}
  

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(1,8,13,"X"))}
fake.f2 <- calc.genoprob(fake.f2, step=5)
out.2dim <- scantwo(fake.f2, method="hk")

# All pairs of chromosomes
summary(out.2dim)

# Chromosome pairs meeting specified criteria
summary(out.2dim, thresholds=c(9.1, 7.1, 6.3, 6.3, 3.3))

# Similar, but ignoring the interaction LOD score in the rule
summary(out.2dim, thresholds=c(9.1, 7.1, Inf, 6.3, 3.3))

# Pairs having largest interaction LOD score, if it's > 4
summary(out.2dim, thresholds=c(0, Inf, 4, Inf, Inf), what="int")

# permutation test to get thresholds; run in two batches
#     and then combined with c.scantwoperm
\dontrun{operm.2dimA <- scantwo(fake.f2, method="hk", n.perm=500)
operm.2dimB <- scantwo(fake.f2, method="hk", n.perm=500)
operm.2dim <- c(operm.2dimA, operm.2dimB)}
\dontshow{strata <- pull.pheno(fake.f2, "sex") + pull.pheno(fake.f2, "pgm")*2
operm.2dim <- scantwo(fake.f2, method="hk", n.perm=100, perm.strata=strata)}

# estimated LOD thresholds
summary(operm.2dim)

# Summary, citing significance levels and so estimating thresholds
#     from the permutation results
summary(out.2dim, perms=operm.2dim, alpha=rep(0.05, 5))

# Similar, but ignoring the interaction LOD score in the rule
summary(out.2dim, perms=operm.2dim, alpha=c(0.05, 0.05, 0, 0.05, 0.05))

# Similar, but also getting genome-scan-adjusted p-values
summary(out.2dim, perms=operm.2dim, alpha=c(0.05, 0.05, 0, 0.05, 0.05),
        pvalues=TRUE)
}

\seealso{ \code{\link[qtl]{scantwo}}, \code{\link[qtl]{plot.scantwo}},
          \code{\link[qtl]{max.scantwo}}, \code{\link[qtl]{condense.scantwo}} }

\keyword{print}
