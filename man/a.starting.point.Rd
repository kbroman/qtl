\name{A starting point}
\docType{package}
\alias{qtl-package}

\title{Introductory comments on R/qtl}
\description{
  A brief introduction to the R/qtl package, with a walk-through of an
  analysis.
}

\section{New to R and/or R/qtl?}{
  \itemize{
    \item In order to use the R/qtl package, you must type (within R)
    \code{library(qtl)}.  You may wish to include this in a
    \code{\link[base:Startup]{.Rprofile}} file.

    \item Documention and several tutorials are available
    at the R archive (\url{https://cran.r-project.org}).

    \item Use the \code{\link[utils]{help.start}} function to start the
    html version of the R help.

    \item Type \code{library(help=qtl)} to get a list of the functions
    in R/qtl.

    \item Use the \code{\link[utils]{example}} function to run examples
    of the various functions in R/qtl.

    \item A tutorial on the use of R/qtl is distributed with the
    package and is also available at \url{http://www.rqtl.org/rqtltour.pdf}.

    \item Download the latest version of R/qtl from the R archive or
    from \url{http://www.rqtl.org}.
  }
}

\section{Walk-through of an analysis}{
Here we briefly describe the use of R/qtl to analyze an experimental
cross.  A more extensive tutorial on its use is distributed with the
package and is also available at \url{http://www.rqtl.org/rqtltour.pdf}.

A difficult first step in the use of most data analysis software is the
import of data.  With R/qtl, one may import data in several different
formats by use of the function \code{\link{read.cross}}.  The
internal data structure used by R/qtl is rather complicated, and is
described in the help file for \code{\link{read.cross}}.  We won't
discuss data import any further here, except to say that the
comma-delimited format (\code{"csv"}) is recommended.  If you have
trouble importing data, send an email to Karl Broman,
\email{kbroman@biostat.wisc.edu}, perhaps attaching examples of your data
files.  (Such data will be kept confidential.)  Also see the sample data
files and code at \url{http://www.rqtl.org/sampledata}.

We consider the example data \code{\link{hyper}}, an experiment on
hypertension in the mouse, kindly provided
by Bev Paigen and Gary Churchill.  Use the \code{\link[utils]{data}}
function to load the data.

\code{data(hyper)}

The \code{\link{hyper}} data set has class \code{"cross"}.  The
function \code{\link{summary.cross}} gives summary information
on the data, and checks the data for internal consistency.  A number
of other utility functions are available; hopefully these are
self-explanatory.

\code{summary(hyper)} \cr
\code{nind(hyper)} \cr
\code{nphe(hyper)} \cr
\code{nchr(hyper)} \cr
\code{nmar(hyper)} \cr
\code{totmar(hyper)}

The function \code{\link{plot.cross}} gives a graphical summary of
the data; it calls \code{\link{plotMissing}} (to plot a matrix
displaying missing genotypes) and \code{\link{plotMap}} (to plot
the genetic maps), and also displays histograms or barplots of the
phenotypes.  The \code{\link{plotMissing}} function can plot
individuals ordered by their phenotypes; you can see that for most
markers, only individuals with extreme phenotypes were genotyped.

\code{plot(hyper)} \cr
\code{plotMissing(hyper)} \cr
\code{plotMissing(hyper, reorder=TRUE)} \cr
\code{plotMap(hyper)}

Note that one marker (on chromosome 14) has no genotype data.  The
function \code{\link{drop.nullmarkers}} removes such markers from
the data.

\code{hyper <- drop.nullmarkers(hyper)} \cr
\code{totmar(hyper)}

The function \code{\link{est.rf}} estimates the recombination
fraction between each pair of markers, and calculates a LOD score for
the test of \eqn{r} = 1/2.  This is useful for identifying markers that
are placed on the wrong chromosome.  Note that since, for these data,
many markers were typed only on recombinant individuals, the pairwise
recombination fractions show rather odd patterns.

\code{hyper <- est.rf(hyper)} \cr
\code{plotRF(hyper)} \cr
\code{plotRF(hyper, chr=c(1,4))}

To re-estimate the genetic map for an experimental cross, use the
function \code{\link{est.map}}.  The function
\code{\link{plotMap}}, in addition to plotting a single map, can
plot the comparison of two genetic maps (as long as they are composed of
the same numbers of chromosomes and markers per chromosome).  The
function \code{\link{replace.map}} map be used to replace the
genetic map in a cross with a new one.

\code{newmap <- est.map(hyper, error.prob=0.01, verbose=TRUE)} \cr
\code{plotMap(hyper, newmap)} \cr
\code{hyper <- replace.map(hyper, newmap)}

The function \code{\link{calc.errorlod}} may be used to assist in
identifying possible genotyping errors; it calculates the error LOD
scores described by Lincoln and Lander (1992).  The
\code{\link{calc.errorlod}} function return a modified version of
the input cross, with error LOD scores included.  The function
\code{\link{top.errorlod}} prints the genotypes with values above a
cutoff (by default, the cutoff is 4.0).

\code{hyper <- calc.errorlod(hyper, error.prob=0.01)} \cr
\code{top.errorlod(hyper)}

The function \code{\link{plotGeno}} may be used to inspect the
observed genotypes for a chromosome, with likely genotyping errors
flagged.

\code{plotGeno(hyper, chr=16, ind=c(24:34, 71:81))}

Before doing QTL analyses, some intermediate calculations need to be
performed.  The function \code{\link{calc.genoprob}} calculates
conditional genotype probabilities given the multipoint marker data.
\code{\link{sim.geno}} simulates sequences of genotypes from their
joint distribution, given the observed marker data.

As with \code{\link{calc.errorlod}}, these functions return a
modified version of the input cross, with the intermediate calculations
included.  The \code{step} argument indicates the density of the grid on
which the calculations will be performed, and determines the density at
which LOD scores will be calculated.

\code{hyper <- calc.genoprob(hyper, step=2.5, error.prob=0.01)} \cr
\code{hyper <- sim.geno(hyper, step=2.5, n.draws=64, error.prob=0.01)}

The function \code{\link{scanone}} performs a genome scan with a
single QTL model.  By default, it performs standard interval mapping
(Lander and Botstein 1989): use of a normal model and the EM algorithm.
If one specifies \code{method="hk"}, Haley-Knott regression is performed
(Haley and Knott 1992).  These two methods require the results from
\code{\link{calc.genoprob}}.

\code{out.em <- scanone(hyper)} \cr
\code{out.hk <- scanone(hyper, method="hk")}

If one specifies \code{method="imp"}, a genome scan is performed by the
multiple imputation method of Sen and Churchill (2001).  This method
requires the results from \code{\link{sim.geno}}.

\code{out.imp <- scanone(hyper, method="imp")}

The output of \code{\link{scanone}} is a data.frame with class
\code{"scanone"}.  The function \code{\link{plot.scanone}} may be
used to plot the results, and may plot up to three sets of results
against each other, as long as they conform appropriately.

\code{plot(out.em)} \cr
\code{plot(out.hk, col="blue", add=TRUE)} \cr
\code{plot(out.imp, col="red", add=TRUE)} \cr
\code{plot(out.hk, out.imp, out.em, chr=c(1,4), lty=1,} \cr
\code{     col=c("blue","red","black"))}

The function \code{\link{summary.scanone}} may be used to list
information on the peak LOD for each chromosome for which the LOD
exceeds a specified threshold.

\code{summary(out.em)} \cr
\code{summary(out.em, threshold=3)} \cr
\code{summary(out.hk, threshold=3)} \cr
\code{summary(out.imp, threshold=3)}

The function \code{\link{max.scanone}} returns the maximum LOD
score, genome-wide.

\code{max(out.em)} \cr
\code{max(out.hk)} \cr
\code{max(out.imp)}

One may also use \code{\link{scanone}} to perform a permutation
test to get a genome-wide LOD significance threshold.

\code{operm.hk <- scanone(hyper, method="hk", n.perm=1000)}

The result has class \code{"scanoneperm"}.  The
\code{\link{summary.scanoneperm}} function may be used to calculate
LOD thresholds.

\code{summary(operm.hk, alpha=0.05)}

The permutation results may also be used in the
\code{\link{summary.scanone}} function to calculate LOD thresholds
and genome-scan-adjusted p-values.

\code{summary(out.hk, perms=operm.hk, alpha=0.05, pvalues=TRUE)}

We should say at this point that the function
\code{\link[base:save]{save.image}} will save your workspace to disk.  You'll
wish you had used this if R crashes.

\code{save.image()}

The function \code{\link{scantwo}} performs a two-dimensional
genome scan with a two-QTL model.  Methods \code{"em"}, \code{"hk"} and
\code{"imp"} are all available.  \code{\link{scantwo}} is
considerably slower than \code{\link{scanone}}, and can require a
great deal of memory.  Thus, you may wish to re-run
\code{\link{calc.genoprob}} and/or \code{\link{sim.geno}} with
a more coarse grid.

\code{hyper <- calc.genoprob(hyper, step=10, err=0.01)} \cr
\code{hyper <- sim.geno(hyper, step=10, n.draws=64, err=0.01)} \cr
\cr
\code{out2.hk <- scantwo(hyper, method="hk")} \cr
\code{out2.em <- scantwo(hyper)} \cr
\code{out2.imp <- scantwo(hyper, method="imp")}

The output is an object with class \code{scantwo}.  The function
\code{\link{plot.scantwo}} may be used to plot the results.  The
upper triangle contains LOD scores for tests of epistasis, while the
lower triangle contains LOD scores for the full model.

\code{plot(out2.hk)} \cr
\code{plot(out2.em)} \cr
\code{plot(out2.imp)}

The function \code{\link{summary.scantwo}} lists the interesting
aspects of the output.  For each pair of chromosomes \eqn{(k,l)}, it
calculates the maximum LOD score for the full model, \eqn{M_f(k,l)}; a
LOD score indicating evidence for a second QTL, allowing for epistasis),
\eqn{M_{fv1}(k,l)}{M_fv1(k,l)}; a LOD score indicating evidence for
epistasis, \eqn{M_i(k,l)}; the LOD score for the additive QTL model,
\eqn{M_a(k,l)}; and a LOD score indicating evidence for a second QTL,
assuming no epistasis, \eqn{M_{av1}(k,l)}{M_av1(k,l)}.

You must provide five LOD thresholds, corresponding to the above five
LOD scores, and in that order.  A chromosome pair is printed if either
(a) \eqn{M_f(k,l) \ge T_f}{M_f(k,l) >= T_f} and (\eqn{M_{fv1}(k,l) \ge
  T_{fv1}}{M_fv1(k,l) >= T_fv1} or \eqn{M_i(k,l) \ge
  T_i}{M_i(k,l) >= T_i}), or (b)
\eqn{M_a(k,l) \ge T_a}{M_a(k,l) >= T_a} and \eqn{M_{av1}(k,l) \ge
  T_{av1}}{M_av1(k,l) >= T_av1}.

\code{summary(out2.em, thresholds=c(6.2, 5.0, 4.6, 4.5, 2.3))} \cr
\code{summary(out2.em, thresholds=c(6.2, 5.0, Inf, 4.5, 2.3))}

In the latter case, the interaction LOD score will be ignored.

The function \code{\link{max.scantwo}} returns the maximum joint
and additive LODs for a two-dimensional genome scan.

\code{max(out2.em)}

Permutation tests may also performed with \code{\link{scantwo}};
it may take a few days of CPU time.  The output is a list containing the
maxima of the above five LOD scores for each of the imputations.

\code{operm2 <- scantwo(hyper, method="hk", n.perm=100)} \cr
\code{summary(operm2, alpha=0.05)}
}

\section{Citing R/qtl}{
  To cite R/qtl in publications, use the Broman et al. (2003) reference
  listed below.
}

\references{
  Broman, K. W. and Sen,
  \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. (2009) \emph{A
  guide to QTL mapping with R/qtl.}  Springer.  \url{http://www.rqtl.org/book}

  Broman, K. W., Wu, H., Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2003) R/qtl: QTL
  mapping in experimental crosses.  \emph{Bioinformatics} \bold{19},
  889--890.

  Haley, C. S. and Knott, S. A. (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  \emph{Heredity} \bold{69}, 315--324.

  Lander, E. S. and Botstein, D. (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps.  \emph{Genetics}
  \bold{121}, 185--199.

  Lincoln, S. E. and Lander, E. S. (1992) Systematic detection of
  errors in genetic linkage data.  \emph{Genomics} \bold{14}, 604--610.

  Sen, \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. and Churchill, G. A. (2001) A statistical framework for quantitative
  trait mapping.  \emph{Genetics} \bold{159}, 371--387.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{models}
