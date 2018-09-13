\name{drop.nullmarkers}
\alias{drop.nullmarkers}

\title{Drop markers without any genotype data}

\description{
  Drop markers, from the data matrices and genetic maps, that have no
  genotype data.
}

\usage{
drop.nullmarkers(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
}

\value{
  The input object, with any markers lacking genotype data removed from
  the genotype data matrices, genetic maps, and, if applicable, any
  derived data (such as produced by \code{\link{calc.genoprob}}).
  (It might be a good idea to re-derive such things after using this
  function.)
}

\examples{
# removes one marker from hyper
data(hyper)
hyper <- drop.nullmarkers(hyper)

# shouldn't do anything to listeria
data(listeria)
listeria <- drop.nullmarkers(listeria)
}

\seealso{ \code{\link{nullmarkers}}, \code{\link{drop.markers}}, \code{\link{clean.cross}},
  \code{\link{geno.table}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{manip}
