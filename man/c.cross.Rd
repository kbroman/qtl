\name{c.cross}
\alias{c.cross}

\title{Combine data for QTL experiments}

\description{
  Concatenate the data for multiple QTL experiments.
}

\usage{
\method{c}{cross}(\dots)
}

\arguments{
 \item{\dots}{A set of objects of class \code{cross}.  See
   \code{\link{read.cross}} for details.  These must all either be of the
   same cross type or be a combination of backcrosses and intercrosses.
   All crosses must have the same number of chromosomes and chromosome
   names, and the same marker orders and positions, though the set of
   markers need not be precisely the same.}
}

\value{
  The concatenated input, as a \code{cross} object. Additional
  columns are added to the phenotype data indicating which cross an
  individual comes from; another column indicates cross type (0=BC,
  1=intercross), if there are crosses of different types.  The crosses
  are not required to have exactly the same set of phenotypes;
  phenotypes with the same names are assumed to be the same.

  If the crosses have different sets of markers, we
  interpolate marker order, but the cM positions of markers that are in
  common between crosses must be precisely the same in the different crosses.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
junk <- fake.f2
junk <- c(fake.f2,junk)
}

\seealso{ \code{\link{subset.cross}} }

\keyword{manip}
