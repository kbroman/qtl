\name{plotPheno}
\alias{plotPheno}

\title{Plot a phenotype distribution}

\description{
  Plots a histogram or barplot of the data for a phenotype from an
  experimental cross.
}

\usage{
plotPheno(x, pheno.col=1, \dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{The phenotype column to plot: a numeric index, or the
   phenotype name as a character string. Alternatively, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
 \item{\dots}{Passed to \code{\link[graphics]{hist}} or
   \code{\link[graphics]{barplot}}.}
}

\value{None.}

\details{
  Numeric phenotypes are displayed as histograms with approximately
  \eqn{2\sqrt{n}}{2 sqrt(n)} bins.  Phenotypes that are factors or that
  have very few unique values are displayed as barplots.
}

\examples{
data(fake.bc)
plotPheno(fake.bc, pheno.col=1)
plotPheno(fake.bc, pheno.col=3)
plotPheno(fake.bc, pheno.col="age")
}

\seealso{ \code{\link{plot.cross}}, \code{\link{plotMap}},
  \code{\link{plotMissing}},
  \code{\link[graphics]{hist}}, \code{\link[graphics]{barplot}} }


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
