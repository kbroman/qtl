\name{plot.MQMall}
\alias{plot.MQMall}

\title{ 
	plot.MQMall - Plotting routine to display the results from a scanMQMall / bootstrapMQM scan 
}
\description{
  Plotting routine to display the results from a scanMQMall scan or a bootstrapMQM, in either a contourmap or a multiple QTL plot from R/QTL plot.scanone
}
\usage{
plot.MQMall(result = NULL, type = "C")
}

\arguments{
  \item{result}{ Result object from scanMQMall }
  \item{type}{ "C" for contourplot, "P" for stacked R/QTL plot, "I" for a heatmapplot  }
}
\details{
  ~~ If necessary, more details than the description above ~~ 
}
\value{
	No return, plotting routine
}
\references{ Danny Arends \email{Danny.Arends@gmail.com} }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link{scanMQMall}} - Function called to do multiple trait analysis
	\item \code{\link{bootstrapMQM}} - Function called to estimate confidence interfalls for MQM
}
}
\examples{
	library(MQMpackage)			#load the MQM library
	data(multitrait)
	result <- scanMQMall(multitrait)
	plot.MQMall(result,"P")
	plot.MQMall(result,"C")
	
}
\references{
\itemize{
	\item Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988), The New S Language. Wadsworth & Brooks/Cole.
}
}
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
\keyword{ hplot }
