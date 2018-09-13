\name{flip.order}
\alias{flip.order}
\title{Flip the orders of markers on a set of chromosomes}

\description{
  Flip the orders of markers on a specified set of chromosome,
  so that the markers will be in the reverse order.
}

\usage{flip.order(cross, chr)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{chr}{Vector indicating the chromosomes to flip. This should be a
 vector of character strings referring to chromosomes by name. A logical
 (TRUE/FALSE) vector may also be used.}
}

\value{
  The input \code{cross} object, but with the marker order on the
  specified chromosomes flipped.
}

\details{
  If the cross contains results from \code{\link{calc.genoprob}},
  \code{\link{sim.geno}}, \code{\link{argmax.geno}}, or
  \code{\link{calc.errorlod}}, those results are also updated.

  Results of \code{\link{est.rf}} and \code{\link{markerlrt}} are
  deleted.
}


\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
fake.f2 <- flip.order(fake.f2, c(1, 5, 13))
}

\seealso{ \code{\link{switch.order}} }

\keyword{manip}
