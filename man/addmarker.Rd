\name{addmarker}
\alias{addmarker}

\title{Add a marker to a cross}

\description{
  Add a marker to a cross object.
}

\usage{
addmarker(cross, genotypes, markername, chr, pos)
}

\arguments{
 \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{genotypes}{Vector of numeric genotypes.}
 \item{markername}{Marker name as character string.}
 \item{chr}{Chromosome ID as character string.}
 \item{pos}{Position of marker, as numeric value.}
}

\details{
  Use this function with caution.  It would be best to incorporate new
  data into a single file to be imported with \code{\link{read.cross}}.

  But if you have genotypes on one or two additional markers that you
  want to add, you might load them with \code{\link[utils]{read.csv}}
  and incorporate them with this function.
}

\value{
  The input \code{cross} object with the single marker added.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\examples{
data(fake.f2)

# genotypes for new marker
gi <- pull.geno(fill.geno(fake.f2))[,"D5M197"]

# add marker to cross
fake.f2 <- addmarker(fake.f2, gi, "D5M197imp", "5", 11)
}

\seealso{ \code{\link{pull.markers}}, \code{\link{drop.markers}} }

\keyword{utilities}
