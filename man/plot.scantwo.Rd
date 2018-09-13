\name{plot.scantwo}
\alias{plot.scantwo}

\title{Plot LOD scores for a two-dimensional genome scan}

\description{
  Plot the results of a two-dimensional, two-QTL genome scan.
}

\usage{
\method{plot}{scantwo}(x, chr, incl.markers=FALSE, zlim, lodcolumn=1,
     lower = c("full", "add", "cond-int", "cond-add", "int"),
     upper = c("int", "cond-add", "cond-int", "add", "full"),
     nodiag=TRUE, contours=FALSE, main, zscale=TRUE, point.at.max=FALSE,
     col.scheme = c("viridis", "redblue","cm","gray","heat","terrain","topo"),
     gamma=0.6, allow.neg=FALSE, alternate.chrid=FALSE, \dots)
}

\arguments{
  \item{x}{An object of class \code{"scantwo"}, as output by
    \code{\link{scantwo}}.}
  \item{chr}{Optional vector indicating the chromosomes to plot.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
  \item{incl.markers}{If FALSE, plot LOD scores on an evenly
    spaced grid (not including the results at the markers).}
  \item{zlim}{A vector of length 2 (optional), indicating the z limits
    for the lower-right and upper-left triangles, respectively.  If one
    number is given, the same limits are used for both triangles.  If
    \code{zlim} is missing, the maximum limits are used for each.}
  \item{lodcolumn}{If the scantwo results contain LOD scores for
    multiple phenotypes, this argument indicates which to use in the
    plot.}
  \item{lower}{Indicates which LOD scores should be plotted in the lower
    triangle.  See the details below.}
  \item{upper}{Indicates which LOD scores should be plotted in the upper
    triangle.  See the details below.}
  \item{nodiag}{If TRUE, suppress the plot of the scanone output
    (which is normally along the diagonal.)}
  \item{contours}{If TRUE, add a contour to the plot at 1.5-LOD below
    its maximum, using a call to \code{\link[graphics]{contour}}.  If a
    numeric vector, contours are drawn at these values below the maximum
    LOD.}
  \item{main}{An optional title for the plot.}
  \item{zscale}{If TRUE, a color scale is plotted at the right.}
  \item{point.at.max}{If TRUE, plot an X at the maximum LOD.}
  \item{col.scheme}{Name of color pallet. The default is "viridis"; see
  Option D at \url{https://bids.github.io/colormap}}
  \item{gamma}{Parameter affecting range of colors when
    \code{col.scheme="gray"} or \code{="redblue"}.}
  \item{allow.neg}{If TRUE, allow the plot of negative LOD scores; in
    this case, the z-limits are symmetric about 0.  This option is
    chiefly to allow a plot of difference between LOD scores from
    different methods, calculated via \code{\link[=-.scantwo]{-.scantwo}}.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome
     axis labels, so that they may be more easily distinguished.}
  \item{\dots}{Ignored at this point.}
}

\value{None.}

\details{
  Uses \code{\link[graphics]{image}} to plot a grid of LOD scores.  The
  particular LOD scores plotted in the upper-left and lower-right
  triangles are selected via \code{upper} and \code{lower},
  respectively.  By default, the upper-left triangle contains the
  epistasis LOD
  scores (\code{"int"}), and the lower-right triangle contains the LOD
  scores for the full model (\code{"full"}).
  The diagonal contains either all zeros or the main effects LOD scores
  (from \code{\link{scanone}}).

  The \code{\link{scantwo}} function calculates, for each pair of
  putative QTLs, \eqn{(q_1,q_2)}{(q1,q2)}, the likelihood under the
  null model \eqn{L_0}{L0}, the likelihood under each of the single-QTL
  models, \eqn{L(q_1)}{L(q1)} and \eqn{L(q_2)}{L(q2)}, the likelihood
  under an additive QTL model, \eqn{L_a(q_1,q_2)}{La(q1,q2)}, and the
  likelihood under a full QTL model (including QTL-QTL interaction),
  \eqn{L_f(q_1,q_2)}{Lf(q1,q2)}.

  The five possible LOD scores that may be plotted are the following.
  The epistasis LOD scores (\code{"int"}) are \eqn{LOD_i = \log_{10}
    L_f(q_1,q_2) - \log_{10} L_a(q_1,q_2)}{LODi = log10 Lf(q1,q2) -
    log10 La(q1,q2)}.

  The full LOD scores (\code{"full"}) are
  \eqn{LOD_f = \log_{10} L_f(q_1,q_2) - \log_{10} L_0}{%
    LODj = log10 Lf(q1,q2) - log10 L0}.

  The additive LOD scores (\code{"add"}) are
  \eqn{LOD_a = \log_{10} L_a(q_1,q_2) - \log_{10} L_0}{%
    LODa = log10 La(q1,q2) - log10 L0}.

  In addition, we may calculate, for each pair of
  chromosomes, the difference between the full LOD score and the
  maximum single-QTL LOD scores for that pair of chromosomes
  (\code{"cond-int"}).

  Finally, we may calculate, for each pair of
  chromosomes, the difference between the additive LOD score and the
  maximum single-QTL LOD scores for that pair of chromosomes
  (\code{"cond-add"}).

  If a color scale is plotted (\code{zscale=TRUE}), the axis on the
  left indicates the scale for the upper-left triangle,
  while the axis on the right indicates the scale for the
  lower-right triangle.  Note that the axis labels can get screwed up
  if you change the size of the figure window; you'll need to redo the
  plot.
}

\section{Output of addpair}{

  \bold{Note} that, for output from \code{\link{addpair}} in which the
  new loci are indicated explicitly in the formula, the summary provided
  by \code{plot.scantwo} is somewhat special.  In particular, the
  \code{lower} and \code{upper} arguments are ignored.

  In the case that the formula used in \code{\link{addpair}} was
  not symmetric in the two new QTL, the x-axis in the plot corresponds
  to the first of the new QTL and the y-axis corresponds to the second
  of the new QTL.
}


\examples{
data(hyper)
\dontshow{hyper <- subset(hyper, chr=c(1,4,6,15))}
hyper <- calc.genoprob(hyper, step=5)
\dontshow{hyper <- calc.genoprob(hyper)}

# 2-d scan by EM and by Haley-Knott regression
out2.em <- scantwo(hyper, method="em")
out2.hk <- scantwo(hyper, method="hk")

# plot epistasis and full LOD scores
plot(out2.em)

# plot cond-int in upper triangle and full in lower triangle
#     for chromosomes 1, 4, 6, 15
plot(out2.em, upper="cond-int", chr=c(1,4,6,15))

# plot cond-add in upper triangle and add in lower triangle
#     for chromosomes 1, 4
plot(out2.em, upper="cond-add", lower="add", chr=c(1,4))

# plot the differences between the LOD scores from Haley-Knott
#     regression and the EM algorithm
plot(out2.hk - out2.em, allow.neg=TRUE)
}

\seealso{ \code{\link{scantwo}},
  \code{\link{summary.scantwo}}, \code{\link{plot.scanone}},
  \code{\link[=-.scantwo]{-.scantwo}} }

\author{Hao Wu; Karl W Broman,
  \email{broman@wisc.edu}; Brian Yandell }

\keyword{hplot}
