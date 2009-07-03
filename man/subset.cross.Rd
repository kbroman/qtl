\name{subset.cross}
\alias{subset.cross}
\alias{[.cross}

\title{Subsetting data for QTL experiment}

\description{
  Pull out a specified set of chromosomes and/or individuals from a
  \code{cross} object.
}

\usage{
\method{subset}{cross}(x, chr, ind, \dots)
\method{[}{cross}(x, chr, ind)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{chr}{Optional vector specifying which chromosomes to keep or discard.
   This may be a logical, numeric, or character string vector.  See
 Details, below.}
 \item{ind}{Optional vector specifying which individuals to keep
   discard. This may be a logical, numeric or chacter string vector.
 See Details, below.} 
 \item{\dots}{Ignored at this point.}
}

\value{The input \code{cross} object, but with only the specified subset
  of the data.}

\details{
  The \code{chr} argument may be a logical vector with length equal to
  the number of chromosomes in the input cross \code{x}.  Alternatively, it
  should be a vector of character strings referring to chromosomes by
  name.  Numeric values are converted to strings.  Refer to chromosomes
  with a preceding \code{-} to have all chromosomes but those
  considered.  

  The \code{ind} argument may be a logical vector with length equal to
  the number of individuals in the input cross \code{x}.  Otherwise, its
  treatment depends on whether the input cross contains individual
  identifiers in the phenotype data (which is a column named either
  \code{id} or \code{ID}).

  If there are no individual identifiers in the input cross, the
  \code{ind} argument should be logical or numeric.  If numeric, the
  values should indicate which individuals to retain (if they are all
  positive) or which to omit (if they are all negative).

  If there are individual identifiers in the input cross and the
  argument \code{ind} is numeric, we first look to see whether all of
  the values match an individual identifier.  If some do not match, we
  treat the vector as we would if the input cross had no individual
  identifiers.  If \code{ind} is numeric and all values match an
  individual identifier, we treat them as such identifiers.

  When the values of \code{ind} are to be treated as individual
  identifiers, we either retain those individuals whose IDs match those
  in \code{ind}, or (in the case that all values in \code{ind} are
  preceded by a \code{-}), we omit those individauls whose IDs match
  those in \code{ind}.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
fake.f2.A <- subset(fake.f2, chr=c("5","13"))
fake.f2.B <- subset(fake.f2, ind = -c(1,5,10))
fake.f2.C <- subset(fake.f2, chr=1:5, ind=1:50)

data(listeria)
y <- pull.pheno(listeria, 1)
listeriaB <- subset(listeria, ind = (!is.na(y) & y < 264))

# individual identifiers
listeria$pheno$ID <- paste("mouse", 1:nind(listeria), sep="")
listeriaC <- subset(listeria, ind=c("mouse1","mouse11","mouse21"))
listeriaD <- subset(listeria, ind=c("-mouse1","-mouse11","-mouse21"))

# you can also use brackets (like matrix with rows=chromosomes and columns=individuals)
temp <- listeria[c("5","13"),]  # chr 5 and 13
temp <- listeria[ , 1:10]       # first ten individuals
temp <- listeria[5, 1:10]       # chr 5 for first ten individuals
}

\seealso{ \code{\link[qtl]{pull.map}}, \code{\link[qtl]{drop.markers}}, \code{\link[qtl]{subset.map}} }

\keyword{manip}
