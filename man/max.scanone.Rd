\name{max.scanone}
\alias{max.scanone}

\title{Maximum peak in genome scan}

\description{
  Print the row of the output from \code{\link{scanone}} that
  corresponds to the maximum LOD, genome-wide.
}

\usage{
\method{max}{scanone}(object, chr, lodcolumn=1, na.rm=TRUE, \dots)
}

\arguments{
 \item{object}{An object of the form output by the function
   \code{\link{scanone}}: a data.frame whose third column is the
   LOD score.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
 \item{lodcolumn}{An integer, indicating which of the LOD score columns
   should be considered in pulling out the peak (these are indexed 1,
   2, \dots).}
 \item{na.rm}{A logical indicating whether missing values should be
   removed.}
 \item{\dots}{Ignored.}
}

\value{
  An object of class \code{summary.scanone}, to be printed by
  \code{print.summary.scanone}.  This is a data.frame with one row,
  corresponding to the maximum LOD peak either genome-wide or for the
  particular chromosome specified.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(listeria)
\dontshow{listeria <- subset(listeria,chr=c(5,13),ind=!is.na(listeria$pheno[,1]))}
listeria <- calc.genoprob(listeria, step=2.5)
out <- scanone(listeria, model="2part", upper=TRUE)
# Maximum peak for LOD(p,mu)
max(out)

# Maximum peak for LOD(p,mu) on chr 5
max(out,chr=5)

# Maximum peak for LOD(p,mu) on chromosomes other than chr 13
max(out,chr="-13")

# Maximum peak for LOD(p)
max(out, lodcolumn=2)

# Maximum peak for LOD(mu)
max(out, lodcolumn=3)
}

\seealso{ \code{\link{scanone}}, \code{\link{plot.scanone}},
          \code{\link{summary.scanone}} }

\keyword{print}
