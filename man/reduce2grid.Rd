\name{reduce2grid}
\alias{reduce2grid}

\title{Reduce to a grid of pseudomarkers.}

\description{
  For high-density marker data, rather than run \code{\link{scanone}} at both the
  markers and at a set of pseudomarkers, we reduce to just
  a set of evenly-spaced pseudomarkers
}

\usage{
reduce2grid(cross)
}

\arguments{
 \item{cross}{An object of class \code{cross}.  See
   \code{\link{read.cross}} for details.}
}

\details{
  Genotype probabilities (from \code{\link{calc.genoprob}}) and/or
  imputations (from \code{\link{sim.geno}}) are subset to a grid of
  pseudomarkers.

  This is so that, in the case of high-density markers, we can do the
  genome scan calculations at a smaller set of points (on an
  evenly-spaced grid, but not at the markers) to save computation time.

  You need to first have run \code{\link{calc.genoprob}} and/or
  \code{\link{sim.geno}}, and you must use \code{stepwidth="fixed"}.

  When plotting results with \code{\link{plot.scanone}}, use
  \code{incl.markers=FALSE}, as the output of \code{\link{scanone}}
  won't include information about the marker locations and so will plot
  tick marks only at the first marker on each chromosome.
}

\value{
  The input \code{cross} object with included genotype probabilities or
  imputations subset to an evenly-spaced grid.
}

\examples{
data(hyper)
hyper <- calc.genoprob(hyper, step=2)
hypersub <- reduce2grid(hyper)

\dontrun{out <- scanone(hypersub)
plot(out, incl.markers=FALSE)}
}

\seealso{ \code{\link{calc.genoprob}},
          \code{\link{sim.geno}}, \code{\link{scanone}}, \code{\link{plot.scanone}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{utilities}
