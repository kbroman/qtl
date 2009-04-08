\name{plot.MQMboot}
\alias{plot.MQMboot}

\title{ plot.MQMboot - Plotting routine to display the results from a bootstrapMQM scan  }

\description{
	Plotting routine to display the results from a bootstrapMQM scan 
}
\usage{
plot.MQMboot(result = NULL,\dots)
}

\arguments{
  \item{result}{ Result object from scanMQMall }
  \item{\dots}{ Extra argument passed to the polyplot methode of R. Brouwer }    
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
	\item \code{\link{bootstrapMQM}} - Function called to estimate confidence interfalls for MQM
}
}
\examples{
#Simulated F2 Population
	library(MQMpackage)
	f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
	data(map10)                                             # Mouse genome
	f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
	f2result <- bootstrapMQM(f2cross)                       # Bootstrap MQM to obtain confidence intervals
	plot.MQMboot(f2result)									# Use the nice bootstrap plotting methodes
}
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
\keyword{ hplot }
