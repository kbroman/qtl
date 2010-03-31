\name{mqmplot.heatmap}
\alias{mqmplot.heatmap}

\title{Plotting routine to display the results from a mqm scan on multiple phenotypes}

\description{
  Plotting routine to display the results from a mqm scan on multiple phenotypes in a heatmap.
}

\usage{
mqmplot.heatmap(cross, results, hidelow=TRUE, directed=TRUE, legend=FALSE)
}

\arguments{
  \item{cross}{
    An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
  \item{results}{ Result object from mqmscanall }
  \item{hidelow}{ Hides the low lodscores (< 3) }
  \item{directed}{ Take direction of QTLs into account (takes more time because of QTL direction calculations }
  \item{legend}{ Add a lengend to the plot ? }
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
    \item \code{\link{mqmpermutation}} - Estimate confidence intervals
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
mqmplot.heatmap(multitrait,result)
}

\keyword{ hplot }
