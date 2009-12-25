\name{mqmloadanalysis}
\alias{mqmloadanalysis}

\title{Enables loading of previously saved MQM analysis from file} 

\description{
  Enables loading of previously saved MQM analysis from file
}

\usage{
mqmloadanalysis(file = "MQM_output.txt")
}

\arguments{
  \item{file}{ File to load MQM analysis from }
}

%\details{
%	Not available yet
%}

\value{
Returns a matrix with QTLscores for each location on the chromosome, which can be plotted using the \code{\link{plot.scanone}} function.
}

\author{ 
  % \dannyauthor
}

\seealso{ 
  \itemize{
    % \input{"inst/docs/mqm/standard_seealso.txt"}
  }
}

\examples{
\dontrun{
qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
result <- mqmscan(cross,file="aaa.txt")			# Do a MQM scan of the genome and save it to aaa.txt
res_reloaded <- mqmloadanalysis("aaa.txt")		#Reads in and plots the file aaa.txt
}
}

\keyword{ IO }
