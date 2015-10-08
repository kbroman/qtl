\name{max.scantwo}
\alias{max.scantwo}

\title{Maximum peak in two-dimensional genome scan}

\description{
  Print the pair of loci with the largest LOD score in the results of
  \code{\link{scantwo}}.
}

\usage{
\method{max}{scantwo}(object, lodcolumn=1,
    what=c("best", "full", "add", "int"),
    na.rm=TRUE, \dots)
}

\arguments{
 \item{object}{An object of class \code{scantwo}, the output of the
   function \code{\link{scantwo}}.}
 \item{lodcolumn}{If the scantwo results contain LOD scores for
   multiple phenotypes, this argument indicates which to use.}
 \item{what}{Indicates for which LOD score the maximum should be reported.}
 \item{na.rm}{Ignored.}
 \item{\dots}{Ignored.}
}

\value{
  An object of class \code{summary.scantwo}, to be printed by
  \code{print.summary.scantwo}, with the pair of positions with the
  maximum LOD score. (Which LOD score is considered is indicated by the
  \code{what} argument.)
}

\details{
  This is very similar to the \code{\link{summary.scantwo}}
  function, though this pulls out one pair of positions.

  If \code{what="best"}, we find the pair of positions at which the LOD
  score for the full model (2 QTL + interaction) is maximized, and then
  also print the positions on that same pair of chromosomes at which the
  additive LOD score is maximized.

  In the other cases, we pull out the pair of positions with the largest
  LOD score; which LOD score is considered is indicated by the
  \code{what} argument.
}

\section{Output of addpair}{

  \bold{Note} that, for output from \code{\link{addpair}} in which the
  new loci are indicated explicitly in the formula, the summary provided
  by \code{max.scantwo} is somewhat special.

  All arguments (except, of course, the input
  \code{object}) are ignored.

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
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
fake.f2 <- calc.genoprob(fake.f2, step=10)
out.2dim <- scantwo(fake.f2, method="hk")
max(out.2dim)
}

\seealso{ \code{\link{scantwo}}, \code{\link{plot.scantwo}},
          \code{\link{summary.scantwo}} }

\keyword{print}
