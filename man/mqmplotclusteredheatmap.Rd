\name{mqmplot.clusteredheatmap}
\alias{mqmplot.clusteredheatmap}

\title{Plot clustered heatmap
of MQM scan on multiple phenotypes}

\description{
  Plot the results from a MQM scan on multiple phenotypes.
}

\usage{
mqmplot.clusteredheatmap(cross, mqmresult, directed=TRUE, legend=FALSE,
                         Colv=NA, scale="none", verbose=FALSE,
                         breaks = c(-100,-10,-3,0,3,10,100),
                         col = c("darkblue","blue","lightblue","yellow",
                                 "orange","red"), \dots)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{mqmresult}{ Result object from mqmscanall, the object needs to be of class \code{mqmmulti} }
  \item{directed}{ Take direction of QTLs into account (takes more time because of QTL direction calculations }
  \item{legend}{ If TRUE, add a legend to the plot }
  \item{Colv}{ Cluster only the Rows, the columns (Markers) should not be clustered }
  \item{scale}{ character indicating if the values should be centered and scaled in either
                the row direction or the column direction, or none. The default "none" }
  \item{verbose}{If TRUE, give verbose output.}
  \item{breaks}{ Color break points for the LOD scores }
  \item{col}{ Colors used between breaks }
  \item{\dots}{Additional arguments passed to \code{\link[stats]{heatmap}}.}
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
cresults <- mqmplot.clusteredheatmap(multitrait,result)
groupclusteredheatmap(multitrait,cresults,10)
}

\keyword{ hplot }
