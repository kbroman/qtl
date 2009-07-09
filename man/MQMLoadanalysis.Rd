\name{MQMLoadanalysis}
\alias{MQMLoadanalysis}

\title{Enables loading of previous MQM analysis from file} 


\description{
  Enables loading of previous MQM analysis from file
}

\usage{
MQMLoadanalysis(file = "MQM_output.txt")
}

\arguments{
  \item{file}{ File to load MQM analysis from }
}

%\details{
%	Not available yet
%}

\value{
Returns a matrix with QTLscores for each location on the chromosome, which can be plotted using the \code{\link[qtl]{plot.scanone}} function.
}

\references{ Danny Arends \email{Danny.Arends@gmail.com} }

\author{ Danny Arends \email{Danny.Arends@gmail.com} }

\seealso{ 
\itemize{
	\item \code{\link[qtl]{scanMQM}} - Function called to do single trait analysis
	\item \code{\link[qtl]{scanall}} - Function called to analytse multiple traits in the cross object
}
}

\examples{
\dontrun{
qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
result <- scanMQM(cross,file="aaa.txt")			# Do a MQM scan of the genome and save it to aaa.txt
res_reloaded <- MQMLoadanalysis("aaa.txt")		#Reads in and plots the file aaa.txt
}
}

\keyword{ IO }
