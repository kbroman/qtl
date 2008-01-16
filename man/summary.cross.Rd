\name{summary.cross}
\alias{summary.cross}
\alias{nind}
\alias{nmar}
\alias{totmar}
\alias{nchr}
\alias{nphe}

\title{Print summary of QTL experiment}

\description{
  Print summary information about a \code{cross} object.
}

\usage{
\method{summary}{cross}(object, \dots)
nind(object)
nmar(object)
totmar(object)
nchr(object)
nphe(object)
}

\arguments{
 \item{object}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{\dots}{Ignored at this point.}
}

\value{
The six functions described here return the following values.
 
  \item{summary.cross}{An object of class \code{summary.cross}
    containing a variety of summary information about the cross (this is
    generally printed automatically).}

  \item{nind}{The number of individuals in the cross.}

  \item{nmar}{A vector containing the number of markers on
  each chromosome.}

  \item{totmar}{The total number of markers.}

  \item{nchr}{The number of chromosomes.}

  \item{nphe}{The number of phenotypes.} 
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
data(fake.f2)
summary(fake.f2)
nind(fake.f2)
nmar(fake.f2)
totmar(fake.f2)
nchr(fake.f2)
nphe(fake.f2)
}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{plot.cross}} }

\keyword{print}
