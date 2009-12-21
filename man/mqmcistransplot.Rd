\name{mqmplot_cistrans}
\alias{mqmplot_cistrans}
\alias{mqmplotcistrans}
\alias{cistransplot}

\title{Plotting routine to show cis and trans QTL information}

\description{
  Plotting routine to show cis/trans QTL information
}

\usage{
	mqmplot_cistrans(x, cross, threshold=5, onlyPEAK=TRUE, highPEAK=FALSE, cisarea=10, pch=22, cex=0.5, \dots)
}

\arguments{
  \item{x}{An object of class \code{mqmmulti}. See \code{\link{mqmall}} for details.} 
  \item{cross}{An object of class \code{cross}. See \code{\link{read.cross}} for details.} 
  \item{threshold}{ Threshold value in LOD, Markers that have a 
	LOD score above this threshold are plotted as small squares 
	(see \code{pch} parameter). The markers with LODscores below
	this threshold are not visible 
  }
  \item{onlyPEAK}{ Plot only the peak markers ? (TRUE/FALSE) 
	(Peak markers are markers that have a QTL likelyhood above
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

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl}}

\note{
If you find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}

\seealso{
	\itemize{
		\item \code{\link{mqmall}} - Function called to analyse multiple traits in the cross object
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

%\references{
%	
%}

\keyword{ utilities }
