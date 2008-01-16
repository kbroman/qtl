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
   \code{\link[qtl]{read.cross}} for details.}
}

\value{
  The input object, with any markers lacking genotype data removed from
  the genotype data matrices, genetic maps, and, if applicable, any
  derived data (such as produced by \code{\link[qtl]{calc.genoprob}}).
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

\seealso{ \code{\link[qtl]{drop.markers}}, \code{\link[qtl]{clean.cross}},
  \code{\link[qtl]{geno.table}} }

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{manip}
