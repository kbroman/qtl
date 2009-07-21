\name{summary.scanone}
\alias{summary.scanone}

\title{Summarize the results of a genome scans}

\description{
  Print the rows of the output from \code{\link[qtl]{scanone}} that
  correspond to the maximum LOD for each chromosome, provided that they
  exceed some specified thresholds.
}

\usage{
\method{summary}{scanone}(object, threshold,
        format=c("onepheno", "allpheno", "allpeaks"),
        perms, alpha, lodcolumn=1, pvalues=FALSE, df=FALSE, \dots)
}

\arguments{
 \item{object}{An object output by the function
   \code{\link[qtl]{scanone}}.}
 \item{threshold}{LOD score thresholds. Only peaks with LOD score above
   this value will be returned. This could be a single number or (for
   formats other than \code{"onepheno"}) a threshold for each LOD score
   column.  If \code{alpha} is specified, \code{threshold} should not be.}
 \item{format}{Format for the output; used only if there is more than
   one LOD score column in the input.  See Details, below.}
 \item{perms}{Optional permutation results used to derive thresholds or
   to calculate genome-scan-adjusted p-values.  This must be consistent
   with the \code{object} input, in that it must have the same number of
   LOD score columns, though it can have just one column of permutation
   results, in which case they are reused for all LOD score columns in
   the \code{\link[qtl]{scanone}} output, \code{object}.
   (These can also be permutation results from
   \code{\link[qtl]{scantwo}}, which permutations for a one-dimensional scan.)}
 \item{alpha}{If perms are included, this is the significance level used
   to calculate thresholds for determining which peaks to pull out.
   If \code{threshold} is specified, \code{alpha} should not be.}
 \item{lodcolumn}{If \code{format="onepheno"}, this indicates the LOD
   score column to focus on.  This should be a single number between 1
   and the number of LOD columns in the object input.}
 \item{pvalues}{If TRUE, include columns with genome-scan-adjusted
   p-values in the results.  This requires that \code{perms} be
   provided.} 
 \item{df}{If TRUE, the degrees of freedom associated with the LOD
   scores are shown.} 
 \item{\dots}{Ignored at this point.}
}

\details{
  This function is used to report loci deemed interesting from a one-QTL
  genome scan (by \code{\link[qtl]{scanone}}).

  For \code{format="onepheno"}, we focus on a single LOD score column,
  indicated by \code{lodcolumn}.  The single largest LOD score peak on
  each chromosome is extracted.  If \code{threshold} is specified, only
  those peaks with LOD meeting the threshold will be
  returned.  If \code{perms} and \code{alpha} are specified, a threshold
  is calculated based on the permutation results in \code{perms} for the
  significance level \code{alpha}.  If neither \code{threshold} nor
  \code{alpha} are specified, the peak on each chromosome is returned.
  Again note that with this format, only the LOD score column indicated
  by \code{lodcolumn} is considered in deciding which chromosomes to
  return, but the LOD scores from other columns, at the position with
  maximum LOD score in the \code{lodcolumn} column, are also returned.

  For \code{format="allpheno"}, we consider all LOD score columns, and
  pull out the position, on each chromosome, showing the largest LOD
  score.  The output thus may contain multiple rows for a chromosome.
  Here \code{threshold} may be a vector of LOD score thresholds, one for
  each LOD score column, in which case only those positions for which a
  LOD score column exceeded its threshold are given.  If
  \code{threshold} is a single number, it is applied to all of the LOD 
  score columns.  If \code{alpha} is specified, it must   be a single
  significance level, applied for all LOD score columns, and again
  \code{perms} must be specified, and these are used to calculate the
  LOD score threshold for the significance level \code{alpha}.

  For \code{format="allpeaks"}, the output will contain, for each
  chromosome, the maximum LOD score for each LOD score column, at the
  position at which it achieved its maximum.  Thus, the output will
  contain no more than one row per chromosome, but will contain the
  position and maximum LOD score for each of the LOD score columns. 
  The arguments \code{threshold} and \code{alpha} may be specified as
  for the \code{"allpheno"} format.  The results for a chromosome are
  returned if at least one of the LOD score columns exceeded its
  threshold. 

  If \code{pvalues=TRUE}, and \code{perms} is specified,
  genome-scan-adjusted p-values are calculated for each LOD score
  column, and there are additional columns in the output containing
  these p-values.

  In the case that X-chromosome specific permutations were performed
  (with \code{perm.Xsp=TRUE} in \code{\link[qtl]{scanone}}), autosome-
  and X-chromosome specific thresholds and p-values are calculated by
  the method in Broman et al. (in press).
}

\value{
  An object of class \code{summary.scanone}, to be printed by
  \code{print.summary.scanone}.  
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=c(2,5,6))}
fake.bc <- calc.genoprob(fake.bc, step=5)

# genome scan by Haley-Knott regression
out <- scanone(fake.bc, method="hk")

# permutation tests
\dontrun{operm <- scanone(fake.bc, method="hk", n.perm=1000)
}\dontshow{operm <- scanone(fake.bc, method="hk", n.perm=100)
}

# peaks for all chromosomes
summary(out)

# results with LOD >= 3
summary(out, threshold=3)

# the same, but also showing the p-values
summary(out, threshold=3, perms=operm, pvalues=TRUE)

# results with LOD meeting the 0.05 threshold from the permutation results
summary(out, perms=operm, alpha=0.05)

# the same, also showing the p-values
summary(out, perms=operm, alpha=0.05, pvalues=TRUE)

##### summary with multiple phenotype results
out2 <- scanone(fake.bc, pheno.col=1:2, method="hk")

# permutations
\dontrun{operm2 <- scanone(fake.bc, pheno.col=1:2, method="hk", n.perm=1000)
}\dontshow{operm2 <- scanone(fake.bc, pheno.col=1:2, method="hk", n.perm=100)
}

# results with LOD >= 2 for the 1st phenotype and >= 1 for the 2nd phenotype
#     using format="allpheno"
summary(out2, thr=c(2, 1), format="allpheno")

# The same with format="allpeaks"
summary(out2, thr=c(2, 1), format="allpeaks")

# The same with p-values
summary(out2, thr=c(2, 1), format="allpeaks", perms=operm2, pvalues=TRUE)

# results with LOD meeting the 0.05 significance level by the permutations
#     using format="allpheno"
summary(out2, format="allpheno", perms=operm2, alpha=0.05)

# The same with p-values
summary(out2, format="allpheno", perms=operm2, alpha=0.05, pvalues=TRUE)

# The same with format="allpeaks"
summary(out2, format="allpeaks", perms=operm2, alpha=0.05, pvalues=TRUE)
}

\seealso{ \code{\link[qtl]{scanone}}, \code{\link[qtl]{plot.scanone}},
  \code{\link[qtl]{max.scanone}}, \code{\link[qtl]{subset.scanone}},
  \code{\link[qtl]{c.scanone}}, \code{\link[qtl]{summary.scanoneperm}}
  \code{\link[qtl]{c.scanoneperm}} }

\references{
  Broman, K. W., Sen, \'S, Owens, S. E., Manichaikul, A.,
  Southard-Smith, E. M. and Churchill G. A.  The X chromosome in
  quantitative trait locus mapping.  \emph{Genetics}, to appear.
}


\keyword{print}
