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
 \item{chr}{A vector specifying which chromosomes to keep or discard.
   This may be a logical vector, a numeric vector, or a vector of
   character strings.} 
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
y <- listeria$pheno[,1]
listeriaB <- subset(listeria, ind = (!is.na(y) & y < 264))
}

\seealso{ \code{\link[qtl]{pull.map}}, \code{\link[qtl]{drop.markers}} }

\keyword{manip}
