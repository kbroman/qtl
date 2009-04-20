\name{plot.nice}
\alias{plot.nice}

\title{ plot.nice - Plotting routine to display the results from a bootstrap or scanall scan  }

\description{
	Plotting routine to display the results from a bootstrap or scanall scan
}
\usage{
plot.nice(result = NULL,\dots)
}

\arguments{
  \item{result}{ Result object from scanall }
  \item{\dots}{ Extra argument passed to the polyplot methode of R. Brouwer. Use legendloc=1 to plot legend outside of plotting area }   
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
	No return, plotting routine
}
\author{ Danny Arends \email{Danny.Arends@gmail.com} ; Rutger Brouwer \email{R.W.W.Brouwer@rug.nl} }
\note{ 	
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link{scanall}} - Function called to analyse multiple traits
	\item \code{\link{bootstrap}} - Function called to estimate confidence interfalls for MQM
}
}
\examples{
#Simulated F2 Population
	f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
	data(map10)                                             # Mouse genome
	f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
	f2result <- bootstrapmqm(f2cross)             # Bootstrap MQM to obtain confidence intervals
	plot.nice(f2result)									# Use the nice plotting methodes
}
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
\keyword{ hplot }
