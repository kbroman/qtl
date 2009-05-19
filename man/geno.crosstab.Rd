\name{geno.crosstab}
\alias{geno.crosstab}

\title{Create table of two-locus genotypes}

\description{
  Create a cross tabulation of the genotypes at a pair of markers.
}

\usage{
geno.crosstab(cross, mname1, mname2, eliminate.zeros=TRUE)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{mname1}{The name of the first marker (as a character
 string). (Alternatively, a vector with the two character strings, in
 which case \code{mname2} should not be given.)}
 \item{mname2}{The name of the second marker (as a character string).}
 \item{eliminate.zeros}{If TRUE, don't show the rows and columns that
 have no data.}
}

\value{
  A matrix containing the number of individuals having each possible
  pair of genotypes.  Genotypes for the first marker are in the rows;
  genotypes for the second marker are in the columns.
}


\examples{
data(hyper)
geno.crosstab(hyper, "D1Mit123", "D1Mit156")
geno.crosstab(hyper, "DXMit22", "DXMit16")
geno.crosstab(hyper, c("DXMit22", "DXMit16"))
}

\seealso{ \code{\link[qtl]{geno.table}},
  \code{\link[qtl]{find.marker}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{utilities}
