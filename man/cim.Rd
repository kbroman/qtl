\name{cim}
\alias{cim}

\title{Composite interval mapping}

\description{
  Composite interval mapping by a scheme from QTL Cartographer:
  forward selection at the markers (here, with filled-in genotype data)
  to a fixed number, followed by interval mapping with the selected
  markers as covariates, dropping marker covariates if they are within
  some fixed window size of the location under test.
}

\usage{
cim(cross, pheno.col=1, n.marcovar=3, window=10,
    method=c("em", "imp", "hk", "ehk"),
    imp.method=c("imp", "argmax"), error.prob=0.0001,
    map.function=c("haldane", "kosambi", "c-v", "morgan"),
    n.perm)
}
\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  One may also give a character string matching
   a phenotype name. Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{n.marcovar}{Number of marker covariates to use.}
 \item{window}{Window size, in cM.}
 \item{method}{Indicates whether to use the EM algorithm,
   imputation, Haley-Knott regression, or the extended Haley-Knott
   method.}
 \item{imp.method}{Method used to impute any missing marker genotype
   data.}
 \item{error.prob}{Genotyping error probability assumed when imputing
   the missing marker genotype data.}
 \item{map.function}{Map function used when imputing the missing marker
   genotype data.}
 \item{n.perm}{If specified, a permutation test is performed rather than
   an analysis of the observed data.  This argument defines the number
   of permutation replicates.}
}


\details{
  We first use \code{\link{fill.geno}} to impute any missing marker
  genotype data, either via a simple random imputation or using the
  Viterbi algorithm.

  We then perform forward selection to a fixed number of markers.  These
  will be used (again, with any missing data filled in) as covariates in
  the subsequent genome scan.
}

\value{
  The function returns an object of the same form as the function
  \code{\link{scanone}}:

  If \code{n.perm} is missing, the function returns the scan results as
  a data.frame with three columns: chromosome, position, LOD score.
  Attributes indicate the names and positions of the chosen marker
  covariates.

  If \code{n.perm} > 0, the function results the results of a
  permutation test: a vector giving the genome-wide maximum LOD score in
  each of the permutations.
}

\references{
  Jansen, R. C. (1993) Interval mapping of multiple quantitative trait
  loci.  \emph{Genetics}, \bold{135}, 205--211.

  Jansen, R. C. and Stam, P. (1994) High resolution of quantitative
  traits into multiple loci via interval mapping.  \emph{Genetics},
  \bold{136}, 1447-1455.

  Zeng, Z. B. (1993) Theoretical basis for separation of multiple linked
  gene effects in mapping quantitative trait loci.
  \emph{Proc. Natl. Acad. Sci. USA}, \bold{90}, 10972--10976.

  Zeng, Z. B. (1994) Precision mapping of quantitative trait loci.
  \emph{Genetics}, \bold{136}, 1457--1468.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(hyper)
hyper <- calc.genoprob(hyper, step=2.5)
\dontshow{hyper <- subset(hyper, chr=c(1,4,6,15))}

out <- scanone(hyper)
out.cim <- cim(hyper, n.marcovar=3)
plot(out, out.cim, chr=c(1,4,6,15), col=c("blue", "red"))

add.cim.covar(out.cim, chr=c(1,4,6,15))
}

\seealso{ \code{\link{add.cim.covar}}, \code{\link{scanone}},
  \code{\link{summary.scanone}}, \code{\link{plot.scanone}},
  \code{\link{fill.geno}}}

\keyword{models}
