\name{summaryScantwoOld}
\alias{summaryScantwoOld}

\title{Summarize the results of a two-dimensional genome scan}

\description{
  Summarize the interesting aspects of the results of
  \code{\link{scantwo}}; this is the version of
  \code{\link{summary.scantwo}} that was included in R/qtl version
  1.03 and earlier.
}

\usage{
summaryScantwoOld(object, thresholds = c(0, 0, 0), lodcolumn=1,
                type = c("joint","interaction"), \dots)
}

\arguments{
 \item{object}{An object of class \code{scantwo}, the output of the
   function \code{\link{scantwo}}.}
 \item{thresholds}{A vector of length three, giving LOD thresholds
   for the joint LOD, interaction LOD and single-QTL conditional
   LOD. Negative threshold values are taken relative to the maximum
   joint, interaction, or individual QTL LOD, respectively.}
  \item{lodcolumn}{If the scantwo results contain LOD scores for
    multiple phenotypes, this argument indicates which to use in the
    summary.}
 \item{type}{Indicates whether to pick peaks with maximal joint or
   interaction LOD.}
 \item{\dots}{Ignored at this point.}
}

\value{
  An object of class \code{summary.scantwo.old}, to be printed by
  \code{print.summary.scantwo.old}.  Pairs of loci meeting
  the specified criteria are printed, with their joint LOD, interaction
  LOD, and the conditional LOD for each locus, along with single-point
  P-values calculated by the \eqn{\chi^2}{chi-square} approximation.
  P-values are printed as \eqn{-\log_{10}(P)}{-log10(P)}. \cr

  If the input \code{scantwo} object does not include the results of
  \code{scanone}, the interaction and conditional LOD thresholds are
  ignored, and all pairs of loci for which the joint LOD exceeds its
  threshold are printed, though without their conditional LOD scores.
}

\details{
  For each pair of chromosomes, the pair of loci for which the
  LOD score (either joint or interaction LOD, according to the argument
  \code{type}) is a maximum is considered.  The pair is printed only if
  its joint LOD score exceeds the joint threshold and either (a) the
  interaction LOD score exceeds its threshold or (b) both of the loci have
  conditional LOD scores that are above the conditional LOD threshold,
  where the conditional LOD score for locus \eqn{q_1}{q1}, \eqn{LOD(q_1
    | q_2)}{LOD(q1 | q2)}, is the \eqn{\log_{10}}{log10} likelihood ratio
  comparing the model with \eqn{q_1}{q1} and \eqn{q_2}{q2} acting
  additively to the model with \eqn{q_2}{q2} alone.

  In the case the results of \code{\link{scanone}} are not
  available, the maximum locus pair for each chromosome is printed
  whenever its joint LOD exceeds the joint LOD threshold.

  The criterion used in this summary is due to Gary Churchill and \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}aunak
  Sen, and deserves careful consideration and possible revision.
}

\author{Hao Wu; Karl W Broman,
  \email{broman@wisc.edu}; Brian Yandell }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(1,6,13))}
fake.f2 <- calc.genoprob(fake.f2, step=5)
out.2dim <- scantwo(fake.f2, method="hk")

# All pairs of loci
summaryScantwoOld(out.2dim)

# Pairs meeting specified criteria
summaryScantwoOld(out.2dim, c(7, 3, 3))

# Pairs with both conditional LODs > 2
summaryScantwoOld(out.2dim,c(0,1000,2))

# Pairs with interaction LOD is above 3
summaryScantwoOld(out.2dim,c(0,3,1000))
}

\seealso{ \code{\link{summary.scantwo}},
  \code{\link{scantwo}}, \code{\link{plot.scantwo}},
          \code{\link{max.scantwo}} }

\keyword{print}
