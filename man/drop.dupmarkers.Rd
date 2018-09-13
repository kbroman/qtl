\name{drop.dupmarkers}
\alias{drop.dupmarkers}

\title{Drop duplicate markers}

\description{
  Drop markers with duplicate names; retaining the first of each set,
  with consensus genotyps
}

\usage{
drop.dupmarkers(cross, verbose=TRUE)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{verbose}{If TRUE, print information on the numbers of genotypes
 and markers omitted.  If > 1, give more detailed information
 on genotypes omitted.}
}

\value{
  The input \code{cross} object, with any duplicate markers omitted
  (except for one).  The marker retained will have consensus genotypes;
  if multiple versions of a marker have different genotypes for an
  individual, they will be replaced by \code{NA}.

  Any derived data (such as produced by \code{\link{calc.genoprob}})
  will be stripped off.
}

\examples{
data(listeria)

listeria <- drop.dupmarkers(listeria)
}

\seealso{ \code{\link{drop.nullmarkers}}, \code{\link{pull.markers}}, \code{\link{drop.markers}},
  \code{\link{summary.cross}}, \code{\link{clean.cross}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{manip}
