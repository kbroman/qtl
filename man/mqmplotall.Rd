\name{mqmplot_multitrait}
\alias{mqmplot_multitrait}

\title{Plotting routine to display the results from a mqm scan on multiple phenotypes}

\description{
  Plotting routine to display the results from a mqm scan on multiple phenotypes. 
  It supports four diffferent visualizations: a contourmap, heatmap, 3D graph or a multiple 
  QTL plot from \code{\link{plot.scanone}}
}

\usage{
mqmplot_multitrait(result, type = "C", theta = 30, phi = 15, \dots)
}

\arguments{
  \item{result}{ Result object from scanall }
  \item{type}{ "C" for contourplot, "P" for stacked R/QTL plot, "I" for a heatmapplot, "D" for a 3D plot  }
  \item{theta}{ Horizontal axis rotation in a 3D plot}
  \item{phi}{ Vertical axis rotation in a 3D plot}
  \item{\dots}{Additional arguments passed to \code{\link[graphics]{plot}}.}
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
    \item \code{\link{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
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
mqmplot_multitrait(result,"P")
mqmplot_multitrait(result,"C")
mqmplot_multitrait(result,"I")
mqmplot_multitrait(result,"D")
}

\keyword{ hplot }
