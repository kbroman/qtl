\name{plotModel}
\alias{plotModel}

\title{Plot a QTL model}

\description{
  Plot a graphical representation of a QTL model, with nodes
  representing QTL and line segments representing pairwise
  interactions.
}

\usage{
plotModel(qtl, formula, circrad.rel=0.25, circrad.abs,
         cex.name=1, chronly=FALSE, order, \dots)
}

\arguments{
  \item{qtl}{A QTL object (as created by \code{\link{makeqtl}}) or
    vector of character strings indicating the names for the QTL.  This
  is also allowed to be a list that contains a component named
  \code{"chr"} (and, optionally, components names \code{"pos"} and
  \code{"formula"}).}
  \item{formula}{Optional formula defining the QTL model.  If missing,
  we look for an attribute \code{"formula"} to the input QTL object or a
  item named \code{"formula"} within the QTL object.}
  \item{circrad.rel}{Radius of the circles that
  indicate the QTL, relative to the distance between the circles.}
  \item{circrad.abs}{Optional radius of the circles that
  indicate the QTL; note that the plotting region will have x- and
  y-axis limits spanning 3 units.}
  \item{cex.name}{Character expansion for the QTL names.}
  \item{chronly}{If TRUE and a formal QTL object is given, only the
  chromosome IDs are used to identify the QTL.}
  \item{order}{Optional vector indicating a permutation of the QTL to
  define where they are to appear in the plot.  QTL are placed around a
  circle, starting at the top and going clockwise.}
  \item{\dots}{Passed to the function \code{plot}.}
}

%\details{
%
%}

\value{None.}

\examples{
# plot a QTL model, using a vector of character strings to define the QTL
plotModel(c("1","4","6","15"), formula=y~Q1+Q2+Q3*Q4)

# plot an additive QTL model
data(hyper)
hyper <- calc.genoprob(hyper)
qtl <- makeqtl(hyper, chr=c(1,4,6,15), pos=c(68.3,30,60,18), what="prob")
plotModel(qtl)

# include an interaction
plotModel(qtl, formula=y~Q1+Q2+Q3*Q4)

# alternatively, include the formula as an attribute to the QTL object
attr(qtl, "formula") <- y~Q1+Q2+Q3*Q4
plotModel(qtl)

# if formula given, the attribute within the object is ignored
plotModel(qtl, y~Q1+Q2+Q3+Q4)

# NULL formula indicates additive QTL model
plotModel(qtl, NULL)

# reorder the QTL in the figure
plotModel(qtl, order=c(1,3,4,2))

# show just the chromosome numbers
plotModel(qtl, chronly=TRUE)
}

\seealso{ \code{\link{stepwiseqtl}}, \code{\link{makeqtl}} }


\author{Karl W Broman, \email{broman@wisc.edu} }

\keyword{hplot}
