\name{groupclusteredheatmap}
\alias{groupclusteredheatmap}

\title{Retrieving groups of traits after clustering}

\description{
  Retrieving groups of traits after clustering, for showing
}

\usage{
  groupclusteredheatmap(cross, clusteredheatmapresult, height)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
  \item{clusteredheatmapresult}{ Result object from \code{\link{mqmplot.clusteredheatmap}} }
  \item{height}{ Height at which to 'cut' the tree, a higher cut-off gives less but larger groups, depending on the tree stucture }
}

\author{ 
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{ 
  \itemize{
    % \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmsetcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
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
