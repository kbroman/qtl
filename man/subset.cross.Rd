\name{subset.cross}
\alias{subset.cross}

\title{Subsetting data for QTL experiment}

\description{
  Pull out a specified set of chromosomes and/or individuals from a
  \code{cross} object.
}

\usage{
\method{subset}{cross}(x, chr, ind, \dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Optional vector specifying which chromosomes to keep.
 This should be a vector of character strings referring to
 chromosomes by name; numeric values are converted to strings.  Refer to
 chromosomes with a preceding \code{-} to have all chromosomes but those
 considered.  A logical (TRUE/FALSE) vector may also be used.} 
 \item{ind}{A vector specifying which individuals to keep or
   discard. This may be a logical or numeric vector.}
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{cross} object, but with only the specified subset
  of the data.}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
fake.f2.A <- subset(fake.f2, chr=c("5","13"))
fake.f2.B <- subset(fake.f2, ind = -c(1,5,10))
fake.f2.C <- subset(fake.f2, chr=1:5, ind=1:50)

data(listeria)
y <- pull.pheno(listeria, 1)
listeriaB <- subset(listeria, ind = (!is.na(y) & y < 264))
}

\seealso{ \code{\link[qtl]{pull.map}}, \code{\link[qtl]{drop.markers}} }

\keyword{manip}
