\name{geno.table}
\alias{geno.table}

\title{Create table of genotype distributions}

\description{
  Create table showing the observed numbers of individuals with each
  genotype at each marker, including P-values from chi-square tests
  for Mendelian segregation.
}

\usage{
geno.table(cross, chr)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{A vector specifying which chromosomes to consider.
   This may be a logical vector, a numeric vector, or a vector of
   character strings.} 
}

\value{
  A matrix containing, for each marker, the number of individuals with
  each possible genotype, as well as the number that were not typed.
  The first column gives the chromosome ID, and the last column gives
  P-values from chi-square tests of Mendelian segregation.  
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
}

\seealso{ \code{\link[qtl]{summary.cross}},
  \code{\link[qtl]{drop.markers}}, \code{\link[qtl]{drop.nullmarkers}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{utilities}
