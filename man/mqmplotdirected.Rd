\name{mqmplot_directedqtl}
\alias{mqmplot_directedqtl}

\title{Plotting routine to display the effects and results from a QTL scan  }

\description{
  Plotting routine to display the effects and results from a mqmscan
}

\usage{
mqmplot_directedqtl(cross, mqmresults, draw = TRUE)
}

\arguments{
  \item{cross}{ 
    % \crossobject
  }
  \item{mqmresults}{ Results from mqmscan of type \code{scanone} }
  \item{draw}{ If TRUE, draw the figure. }  
}

\value{
  Returns a scanone object, with added the effectsign calculated 
  internally by the function effect.scan. For more info on the 
  scanone object see: \code{\link{scanone}}
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
mqmplot_directedqtl(f2cross,f2result)
}

\keyword{ hplot }
