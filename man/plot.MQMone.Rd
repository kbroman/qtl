\name{plot.MQMone}
\alias{plot.MQMone}

\title{ plot.MQMone - Plotting routine to display the results from a MQMscan  }

\description{
  Plotting routine to display the results from a MQMscan
}
\usage{
plot.MQMone(result = NULL, result2=NULL, extended = 0,\dots)
}

\arguments{
  \item{result}{ Results from scanMQM of type scanone }
  \item{result2}{ Results from scanMQM of type scanone }  
  \item{extended}{ Extended plotting of the information content }
  \item{\dots}{ Extra argument passed to the plot.scanone }  
}

\details{
  ~~ If necessary, more details than the description above ~~
}

\value{
	No return, plotting routine
}

\author{ Danny Arends \email{Danny.Arends@gmail.com}  }

\note{ 	
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
  \itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
  }
}
\examples{
#Simulated F2 Population
	library(MQMpackage)
	f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
	data(map10)                                             # Mouse genome
	f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
	f2result <- scanMQM(f2cross)                    	   # Do a MQM scan of the genome
	plot.MQMone(f2result)									#Use our fancy plotting routine
}

\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
\keyword{ hplot }
