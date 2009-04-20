\name{plot.all}
\alias{plot.all}

\title{ 
	plot.all - Plotting routine to display the results from a scanall / bootstrap scan 
}
\description{
  Plotting routine to display the results from a scanall scan or a bootstrap, in either a contourmap or a multiple QTL plot from R/QTL plot.scanone
}
\usage{
plot.all(result = NULL, type = "C", theta = 30, phi = 15)
}

\arguments{
  \item{result}{ Result object from scanall }
  \item{type}{ "C" for contourplot, "P" for stacked R/QTL plot, "I" for a heatmapplot, "D" for a 3D plot  }
  \item{theta}{ Horizontal axis rotation in a 3D plot}
  \item{phi}{ Vertical axis rotation in a 3D plot}
}
\details{

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
	\item \code{\link{scanall}} - Function called to analyse multiple traits in the cross object
	\item \code{\link{bootstrap}} - Function called to estimate confidence interfalls for MQM
}
}
\examples{
	data(multitrait)
	result <- mqmall(multitrait)
	plot.all(result,"P")
	plot.all(result,"C")
	
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
