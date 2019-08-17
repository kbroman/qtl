\name{mqmplot.heatmap}
\alias{mqmplot.heatmap}

\title{Heatmap of a genome of MQM scan on multiple phenotypes}

\description{
  Plotting routine to display a heatmap of results obtained from a multiple-QTL model on multiple phenotypes (the
  output of \code{\link{mqmscanall}})
}

\usage{
mqmplot.heatmap(cross, result, directed=TRUE, legend=FALSE, breaks =
c(-100,-10,-3,0,3,10,100), col =
c("darkblue","blue","lightblue","yellow","orange","red"), \dots)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{result}{ Result object from mqmscanall, the object needs to be of class \code{mqmmulti} }
  \item{directed}{ Take direction of QTLs into account (takes more time because of QTL direction calculations }
  \item{legend}{ If TRUE, add a legend to the plot }
  \item{breaks}{ Color break points for the LOD scores }
  \item{col}{ Colors used between breaks }
  \item{\dots}{Additional arguments passed to the \code{\link{image}} function }
}

\author{
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{
  \itemize{
    % \input{"inst/doc/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{https://rqtl.org/tutorials/MQM-tour.pdf}
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
mqmplot.heatmap(multitrait,result)
}

\keyword{ hplot }
