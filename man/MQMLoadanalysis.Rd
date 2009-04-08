\name{MQMLoadanalysis}
\alias{MQMLoadanalysis}

\title{ 
	MQMLoadanalysis - Enables loading of previous MQM analysis from file 
}
\description{
  Enables loading of previous MQM analysis from file
}
\usage{
MQMLoadanalysis(file = "MQM_output.txt")
}

\arguments{
  \item{file}{ File to load MQM analysis from }
}
\details{
	Not available yet
}
\value{
	Returns a matrix with QTLscores for each location on the chromosome, which can be plotted using the \link[qtl]{plot.scanone} function of R/QTL.
}
\references{ Danny Arends \email{Danny.Arends@gmail.com} }
\author{ Danny Arends \email{Danny.Arends@gmail.com} }
\note{ 	
	If u find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}
\seealso{ 
\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
}
\examples{
library(MQMpackage)
qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
result <- scanMQM(cross,file="aaa.txt")			# Do a MQM scan of the genome and save it to aaa.txt
res_reloaded <- MQMLoadanalysis("aaa.txt")		#Reads in and plots the file aaa.txt
}

\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
