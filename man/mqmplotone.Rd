\name{mqmplot_one}
\alias{mqmplot_one}

\title{Plotting routine to display the results from a QTL scan  }

\description{
  Plotting routine to display the results from a mqmscan
}

\usage{
mqmplot_one(result, extended = 0 ,\dots)
}

\arguments{
  \item{result}{ Results from mqmscan of type \code{scanone}} 
  \item{extended}{ Extended plotting of the information content }
  \item{\dots}{ Extra argument passed to the \code{\link{plot.scanone}} }  
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
#Simulated F2 Population
f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
data(map10)                                             # Mouse genome
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
f2cross <- mqmaugment(f2cross)
f2result <- mqmscan(f2cross)                    	# Do a MQM scan of the genome
mqmplot_one(f2result)					# Use our fancy plotting routine
}

\keyword{ hplot }
