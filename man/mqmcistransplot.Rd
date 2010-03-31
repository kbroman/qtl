\name{mqmplot.cistrans}
\alias{mqmplot.cistrans}

\title{Plotting routine to show cis and trans QTL information}

\description{
  Plotting routine to show cis/trans QTL information
}

\usage{
	mqmplot.cistrans(x, cross, threshold=5, onlyPEAK=TRUE, highPEAK=FALSE, cisarea=10, pch=22, cex=0.5, \dots)
}

\arguments{
  \item{x}{An object of class \code{mqmmulti}. See \code{\link{mqmscanall}} for details.} 
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
  \item{threshold}{ Threshold value in LOD, Markers that have a 
	LOD score above this threshold are plotted as small squares 
	(see \code{pch} parameter). The markers with LODscores below
	this threshold are not visible 
  }
  \item{onlyPEAK}{ Plot only the peak markers ? (TRUE/FALSE) 
	(Peak markers are markers that have a QTL likelihood above
	\code{threshold} and higher than other markers in the same section) 
  }
  \item{highPEAK}{ Highlight peak markers ? (TRUE/FALSE). 
	When using this option peak markers (the marker with the 
	highest LOD score in a region above the threshold gets 
	an 25\% increase in size and is displayed in red) 
  }
  \item{cisarea}{ Adjust the two green lines around the line y=x }
  \item{pch}{ What kind of character is used in plotting of the figure (Default: 22, small square) }
  \item{cex}{ Size of the points plotted (default to 0.5 half of the original size)}
  \item{\dots}{ Extra paramters will be passed to points }
}

%\details{
%
%}

\value{
	Plotting routine, so no return
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
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
	}
}

\examples{
	data(multitrait)
\dontshow{multitrait <- subset(multitrait, chr=1, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))}
	data(locations)
	multiloc <- addloctocross(multitrait,locations)
    multiloc <- calc.genoprob(multiloc)
	results <- scanall(multiloc, method="hk")
	mqmplot_cistrans(results, multiloc, 5, FALSE, TRUE)
}

\keyword{ hplot }
