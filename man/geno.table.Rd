\name{geno.table}
\alias{geno.table}

\title{Create table of genotype distributions}

\description{
  Create table showing the observed numbers of individuals with each
  genotype at each marker, including P-values from chi-square tests
  for Mendelian segregation.
}

\usage{
geno.table(cross, chr, scanone.output=FALSE)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{scanone.output}{If TRUE, give result in the form output by
 \code{\link{scanone}}, so that one may use \code{\link{plot.scanone}}, etc.}
}

\value{
  If \code{scanone.output=FALSE}, the output is a matrix containing, for
  each marker, the number of individuals with each possible genotype, as
  well as the number that were not typed.  The first column gives the
  chromosome ID, and the last column gives P-values from chi-square
  tests of Mendelian segregation.

  If \code{scanone.output=TRUE}, the output is of the form produced by
  \code{\link{scanone}}, with the first two columns being chromosome IDs
  and cM positions of the markers.  The third column is
  \eqn{-\log_{10}(P)}{-log10P} from chi-square tests of Mendelian
  segregation.  The fourth column is the proportion of missing data.
  The remaining columns are the proportions of the different genotypes
  (among typed individuals).
}

\details{
  The P-values are obtained from chi-square tests of Mendelian
  segregation.  In the case of the X chromosome, the sexes and cross
  directions are tested separately, and the chi-square statistics
  combined, and so the test is of whether any of the groups show
  deviation from Mendel's rules.
}

\examples{
data(listeria)
geno.table(listeria)

geno.table(listeria, chr=13)

gt <- geno.table(listeria)
gt[gt$P.value < 0.01,]

out <- geno.table(listeria, scanone.output=TRUE)
plot(out)
plot(out, lod=2)
}

\seealso{ \code{\link{summary.cross}},
  \code{\link{drop.markers}}, \code{\link{drop.nullmarkers}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{utilities}
