\name{mqmplot.multitrait}
\alias{mqmplot.multitrait}

\title{Plot the results from a genomescan using a multiple-QTL model on multiple phenotypes}

\description{
  Plotting routine to display the results from a  multiple-QTL model on
  multiple phenotypes.  It supports four different visualizations: a
  contourmap, heatmap, 3D graph or a multiple QTL plot created by using
  \code{\link{plot.scanone}} on the \code{mqmmulti} object
}

\usage{
mqmplot.multitrait(result, type=c("lines","image","contour","3Dplot"),
                   group=NULL, meanprofile=c("none","mean","median"),
                   theta=30, phi=15, \dots)
}

\arguments{
  \item{result}{ Result object from \code{\link{mqmscanall}} }
  \item{type}{ Selection of the plot method to visualize the data: "lines" (defaut plotting option), "image", "contour" and "3Dplot" }
  \item{group}{ A numeric vector indicating which traits to plot. NULL means no grouping }
  \item{meanprofile}{ Plot a mean/median profile from the group selected }
  \item{theta}{ Horizontal axis rotation in a 3D plot }
  \item{phi}{ Vertical axis rotation in a 3D plot }
  \item{\dots}{Additional arguments passed to \code{\link[graphics]{plot}}.}
}

\author{
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{
  \itemize{
    % \input{"inst/doc/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{http://rqtl.org/tutorials/MQM-tour.pdf}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Set cofactors using marker density
    \item \code{\link{mqmsetcofactors}} - Set cofactors at fixed locations
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/doc/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
data(multitrait)
\dontshow{multitrait <- subset(multitrait, chr=1:2, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
multitrait$pheno <- multitrait$pheno[,1:3]}
multitrait <- fill.geno(multitrait) # impute missing genotype data
result <- mqmscanall(multitrait, logtransform=TRUE)
mqmplot.multitrait(result,"lines")
mqmplot.multitrait(result,"contour")
mqmplot.multitrait(result,"image")
mqmplot.multitrait(result,"3Dplot")
}

\keyword{ hplot }
