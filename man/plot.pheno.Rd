\name{plot.pheno}
\alias{plot.pheno}

\title{Plot a phenotype distribution}

\description{
  Plots a histogram or barplot of the data for a phenotype from an
  experimental cross. 
}

\usage{
\method{plot}{pheno}(x, pheno.col=1, \dots)
}

\arguments{
 \item{x}{An object of class \code{cross}.  See
   \code{\link[qtl]{read.cross}} for details.}
 \item{pheno.col}{The phenotype column to plot: a numeric index, or the
   phenotype name as a character string.}
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
plot(fake.bc, pheno.col=1)
plot(fake.bc, pheno.col=3)
}

\seealso{ \code{\link[qtl]{plot.cross}}, \code{\link[qtl]{plot.map}}, 
  \code{\link[qtl]{plot.missing}}, 
  \code{\link[graphics]{hist}}, \code{\link[graphics]{barplot}} }


\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\keyword{hplot}
