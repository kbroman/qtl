\name{plotLodProfile}
\alias{plotLodProfile}

\title{Plot 1-d LOD profiles for a multiple QTL model}

\description{
  Use the results of \code{\link{refineqtl}} to plot
  one-dimensional LOD profiles for each QTL.
}

\usage{
plotLodProfile(qtl, chr, incl.markers=TRUE, gap=25, lwd=2, lty=1, col="black",
               qtl.labels=TRUE, mtick=c("line", "triangle"),
               show.marker.names=FALSE, alternate.chrid=FALSE,
               add=FALSE, showallchr=FALSE, labelsep=5, \dots)

}			    
\arguments{
  \item{qtl}{An object of class \code{"qtl"}; must have been produced by
    \code{\link{refineqtl}} using \code{keeplodprofiles=TRUE}.}
  \item{chr}{Optional vector indicating the chromosomes to plot.
  This should be a vector of character strings referring to chromosomes
  by name; numeric values are converted to strings.  Refer to
  chromosomes with a preceding \code{-} to have all chromosomes but
  those considered.  A logical (TRUE/FALSE) vector may also be used.}
  \item{incl.markers}{Indicate whether to plot line segments at the
  marker locations.}
  \item{gap}{Gap separating chromosomes (in cM).}
  \item{lwd}{Line widths for each QTL trace (length 1 or the number of QTL).}
  \item{lty}{Line types for each QTL trace (length 1 or the number of QTL).}
  \item{col}{Line col for each QTL trace (length 1 or the number of QTL).}
  \item{qtl.labels}{If TRUE, place a label on each QTL trace.}
  \item{mtick}{Tick mark type for markers (line segments or
    upward-pointing triangels).}
  \item{show.marker.names}{If TRUE, show the marker names along the x axis.}
  \item{alternate.chrid}{If TRUE and more than one chromosome is
    plotted, alternate the placement of chromosome 
     axis labels, so that they may be more easily distinguished.}
  \item{add}{If TRUE, add curves to a current plot.}
  \item{showallchr}{If FALSE (the default), only show the chr with a QTL}
  \item{labelsep}{If \code{qtl.labels=TRUE}, separation between peak LOD
    and QTL label, as percent of the height of the plot.}
  \item{\dots}{Passed to the function \code{\link[graphics]{plot}} when it
    is called.}
}
\value{None.}

\details{
  The function plots LOD profiles in the context of a multiple QTL
  model, using a scheme best described in Zeng et al. (2000).  The
  position of each QTL is varied, keeping all other loci fixed.  If a
  QTL is isolated on a chromosome, the entire chromosome is scanned; if
  there are additional linked QTL, the position of a QTL is scanned over
  the largest interval possible without allowing the order of QTLs along
  a chromosome to change.  At each position for the QTL being scanned,
  we calculate a LOD score comparing the full model, with the QTL of
  interest at that particular position (and all others at their fixed
  positions) to the model with the QTL of interest (and any
  interactions that include that QTL) omitted.  

  Care should be take regarding the arguments \code{lwd}, \code{lty},
  and \code{col}; if vectors are given, they should be in the order of
  the QTL within the object, which may be different than the order in
  which they are plotted.  (The LOD profiles are sorted by chromosome
  and position.)
}

\examples{
data(fake.bc)
\dontshow{fake.bc <- subset(fake.bc, chr=c(2,5))}
fake.bc <- calc.genoprob(fake.bc, step=2)
qtl <- makeqtl(fake.bc, chr=c(2,5), pos=c(32.5, 17.5), what="prob")

out <- scanone(fake.bc, method="hk")

# refine QTL positions and keep LOD profiles
rqtl <- refineqtl(fake.bc, qtl=qtl, method="hk", keeplodprofile=TRUE)

# plot the LOD profiles
plotLodProfile(rqtl)

# add the initial scan results, for comparison
plot(out, add=TRUE, chr=c(2,5), col="red")
}

\seealso{ \code{\link{refineqtl}},
  \code{\link{makeqtl}}, \code{\link{scanqtl}} }

\author{Karl W Broman, \email{broman@wisc.edu} }

\references{
  Zeng  Z.-B., Liu, J., Stam, L. F., Kao, C.-H., Mercer, J. M. and
  Laurie, C. C. (2000) Genetic architecture of a morphological shape
  difference between two Drosophila species. \emph{Genetics} \bold{154},
  299--310. 
}

\keyword{hplot}
