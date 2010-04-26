\name{find.pheno}
\alias{find.pheno}

\title{Find column number for a particular phenotype}

\description{
  Find the column number corresponding to a particular
  phenotype name.
}

\usage{find.pheno(cross, pheno)}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{pheno}{Vector of phenotype names (as character strings).}
}

\value{
  A vector of numbers, corresponding to the column numbers of the
  phenotype in the input cross with the specified names.
}

%\details{
%}

\author{Brian Yandell}

\examples{
data(fake.bc)
find.pheno(fake.bc, "sex")
}

%\seealso{ }

\keyword{utilities}
