\name{CisTransPlot}
\alias{CisTransPlot}

\title{Plotting routine to show cis and trans QTL information}

\description{
  Plotting routine to show cis/trans QTL information
}

\usage{
	CisTransPlot(x, cross, threshold=5, onlyPEAK=TRUE, highPEAK=FALSE, cisarea=10, pch=22, cex=0.5, \dots)
}

\arguments{
  \item{x}{An object of class \code{mqmmulti}. See \code{\link[qtl]{scanall}} for details.} 
  \item{cross}{An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details.} 
  \item{threshold}{ Threshold value in LOD, QTLscores below this threshold are ignored }
  \item{onlyPEAK}{ Plot only the peak markers ? (TRUE/FALSE) (Peak markers are markers that have a QTL likelyhood above \code{threshold} and higher than other markers in the same section) }
  \item{highPEAK}{ highlight peak markers ? (TRUE/FALSE) (peak markers get 25\% increase in size and are colored red) }
  \item{cisarea}{ Adjust the two green lines around the y=x }
  \item{pch}{ What kind of point is plotted (22=small square) }
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
		\item \code{\link[qtl]{addloctocross}} - Function called to analytse multiple traits in the cross object
	}
}

\examples{
	data(multitrait)
\dontshow{multitrait <- subset(multitrait, chr=1, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))}
	data(locations)
	multiloc <- addloctocross(multitrait,locations)
        multiloc <- calc.genoprob(multiloc)
	results <- scanall(multiloc, method="hk")
	CisTransPlot(results, multiloc, 5, FALSE, TRUE)
}

%\references{
%	
%}

\keyword{ utilities }
