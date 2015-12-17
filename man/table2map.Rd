\name{table2map}
\alias{table2map}

\title{Convert a table of marker positions to a map object.}

\description{
  Convert a data frame with marker positions to a map object.
}

\usage{table2map(tab)}

\arguments{
 \item{tab}{A data frame with two columns: chromosome and position.
  The row names are the marker names.}
}

\value{
 A \code{map} object: a list whose components (corresponding
 to chromosomes) are vectors of marker positions.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\examples{
tab <- data.frame(chr=c(1,1,1,1,2,2,2,2,3,3,3,3),
                  pos=c(0,2,4,8,0,2,4,8,0,2,4,8))
rownames(tab) <- paste0("marker", 1:nrow(tab))

map <- table2map(tab)
}

\seealso{ \code{\link{map2table}}, \code{\link{pull.map}}, \code{\link{est.map}} }


\keyword{utilities}
