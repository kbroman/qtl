\name{plotPXG}
\alias{plotPXG}

\title{Plot phenotypes versus marker genotypes}

\description{
  Plot the phenotype values versus the genotypes at a marker or markers.
}

\usage{
plotPXG(x, marker, pheno.col=1, jitter=1, infer=TRUE,
         pch, ylab, main, col, \dots)
}

\arguments{
  \item{x}{An object of class \code{cross}.  See
    \code{\link{read.cross}} for details. }
  \item{marker}{Marker name (a character string; can be a vector).}
  \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  One may also give a character string matching
   a phenotype name. Finally, one may give a numeric vector of
   phenotypes, in which case it must have the length equal to the number
   of individuals in the cross, and there must be either non-integers or
   values < 1 or > no. phenotypes; this last case may be useful for studying
   transformations.}
  \item{jitter}{A positive number indicating how much to spread out the
    points horizontally.  (Larger numbers correspond to greater spread.)}
  \item{infer}{If TRUE, missing genotypes are filled in with a single
    random imputation and plotted in red; if FALSE, only individuals
    typed at the specified marker are plotted.}
  \item{pch}{Plot symbol.}
  \item{ylab}{Label for y-axis.}
  \item{main}{Main title for the plot.  If missing, the names of the
    markers are used.}
  \item{col}{A vector of colors to use for the confidence intervals
    (optional).}
  \item{\dots}{Passed to \code{plot}.}
}

\details{
  Plots the phenotype data against the genotypes at the specified
  marker.  If infer=TRUE, the genotypes of individuals that were not
  typed is inferred based the genotypes at linked markers via a single
  imputation from \code{\link{sim.geno}}; these points are plotted
  in red.  For each genotype, the phenotypic mean is plotted, with error
  bars at \eqn{\pm}{+/-} 1 SE.
}

\value{
A data.frame with initial columns the marker genotypes, then the
phenotype data, then a column indicating whether any of the marker
genotypes were inferred (1=at least one genotype inferred, 0=none were
inferred).
}

\examples{
data(listeria)
mname <- find.marker(listeria, 5, 28) # marker D5M357
plotPXG(listeria, mname)

mname2 <- find.marker(listeria, 13, 26) # marker D13Mit147
plotPXG(listeria, c(mname, mname2))
plotPXG(listeria, c(mname2, mname))

# output of the function contains the raw data
output <- plotPXG(listeria, mname)
head(output)

# another example
data(fake.f2)
mname <- find.marker(fake.f2, 1, 37) # marker D1M437
plotPXG(fake.f2, mname)

mname2 <- find.marker(fake.f2, "X", 14) # marker DXM66
plotPXG(fake.f2, mname2)

plotPXG(fake.f2, c(mname,mname2))
plotPXG(fake.f2, c(mname2,mname))

}

\seealso{ \code{\link{find.marker}}, \code{\link{effectplot}},
  \code{\link{find.flanking}}, \code{\link{effectscan}} }


\author{Karl W Broman, \email{broman@wisc.edu}; Brian Yandell }

\keyword{hplot}
