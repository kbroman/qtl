\name{map2table}
\alias{map2table}

\title{Convert genetic map from list to table.}

\description{
  Convert a map object (as a list) to a table (as a data frame).
}

\usage{map2table(map, chr)}

\arguments{
 \item{map}{A \code{map} object: a list whose components (corresponding
 to chromosomes) are either vectors of marker positions or matrices with
 two rows of sex-specific marker positions.}
 \item{chr}{Optional vector indicating the chromosomes to consider.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
}

\value{
  A data frame with two or three columns: chromosome and sex-averaged
  position, or chromosome, female position, and male position.

  The row names are the marker names.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
map <- pull.map(fake.f2)
map_as_tab <- map2table(map)
}

\seealso{ \code{\link{table2map}}, \code{\link{pull.map}}, \code{\link{est.map}} }


\keyword{utilities}
