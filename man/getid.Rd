\name{getid}
\alias{getid}

\title{Pull out the individual identifiers from a cross}

\description{
  Pull out the individual identifiers from a cross object.
}

\usage{getid(cross)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
}

\value{
  A vector of individual identifiers, pulled from the phenotype data (a
  column named \code{id} or \code{ID}).

  If there are no such identifiers in the cross, the function returns
  \code{NULL}.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)

# create an ID column
fake.f2$pheno$id <- paste("ind", sample(nind(fake.f2)), sep="")

getid(fake.f2)
}

\seealso{ \code{\link{subset.cross}}, \code{\link{top.errorlod}} }


\keyword{utilities}
