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
   \code{\link{read.cross}} for details.}
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

  If the \code{ind} argument is a logical vector
  (\code{TRUE}/\code{FALSE}), it should have length equal to the number
  of individuals in the input cross \code{x}.  The individuals with
  corresponding \code{TRUE} values are retained.

  If the \code{ind} argument is numeric, it should have values either
  between 1 and the number of individuals in the input cross \code{x}
  (in which case these individuals will be retained),
  or it should have values between \code{-1} and \code{-n}, where
  \code{n} is the number of individuals in the input cross \code{x}, in
  which case all \emph{except} these individuals will be retained.

  If the input cross object \code{x} contains individual identifiers (a
  phenotype column labeled \code{"id"} or \code{"ID"}), and if the
  \code{ind} argument contains character strings, then these will be
  matched against the individual identifiers.
  If all values in \code{ind} are
  preceded by a \code{-}), we omit those individuals whose IDs match
  those in \code{ind}.  Otherwise, we retain those individuals whose IDs
  match those in \code{ind}.
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

\seealso{ \code{\link{pull.map}}, \code{\link{drop.markers}}, \code{\link{subset.map}} }

\keyword{manip}
