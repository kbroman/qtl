\name{comparecrosses}
\alias{comparecrosses}
\title{Compare two cross objects}

\description{
  Verify that two objects of class \code{cross} have identical classes,
  chromosomes, markers, genotypes, genetic maps, and phenotypes.
}

\usage{comparecrosses(cross1, cross2, tol=1e-5)}

\arguments{
 \item{cross1}{An object of class \code{cross} (must be an
   intercross). See \code{\link{read.cross}} for details.}
 \item{cross2}{An object of class \code{cross} (must be an
   intercross). See \code{\link{read.cross}} for details.}
 \item{tol}{Tolerance value for comparing genetic map positions and
   numeric phenotypes.}
}

\value{None.}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(listeria)
comparecrosses(listeria, listeria)
}

\seealso{ \code{\link{summary.cross}} }

\keyword{utilities}
