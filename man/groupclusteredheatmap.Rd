\name{groupclusteredheatmap}
\alias{groupclusteredheatmap}

\title{Retrieving groups of traits after clustering}

\description{
  Retrieving groups of clustered traits from the output of mqmplot.clusteredheatmap.
}

\usage{
  groupclusteredheatmap(cross, clusteredheatmapresult, height)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{clusteredheatmapresult}{ Resultint dendrogram object from \code{\link{mqmplot.clusteredheatmap}} }
  \item{height}{ Height at which to 'cut' the dendrogram, a higher cut-off gives less but larger groups.
  Height represents the maximum distance between two traits clustered together using hclust. the 'normal'
  behaviour of bigger groups when using a higher heigh cut-off depends on the tree stucture and the amount
  of traits clustered using \code{\link{mqmplot.clusteredheatmap}} }
}

\value{
  A list containing groups of traits which were clustered together with a distance less that \code{height}
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
cresults <- mqmplot.clusteredheatmap(multitrait,result)
groupclusteredheatmap(multitrait,cresults,10)
}

\keyword{utilities}
