\name{plotMQMone}
\alias{plotMQMone}

\title{Plotting routine to display the results from a QTL scan  }

\description{
  Plotting routine to display the results from a MQMscan
}

\usage{
plotMQMone(result, result2, extended = 0 ,\dots)
}

\arguments{
  \item{result}{ Results from scanMQM of type \code{scanone}}
  \item{result2}{ Results from scanMQM of type \code{scanone} }  
  \item{extended}{ Extended plotting of the information content }
  \item{\dots}{ Extra argument passed to the \code{\link[qtl]{plot.scanone}} }  
}

%\details{
%  ~~ If necessary, more details than the description above ~~
%}

\value{None.}


\author{ Danny Arends \email{Danny.Arends@gmail.com}  }

\seealso{ 
  \itemize{
	\item \code{\link[qtl]{scanMQM}} - Function called to do single trait analysis
  }
}

\examples{
#Simulated F2 Population
f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
data(map10)                                             # Mouse genome
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
f2cross <- MQMaugment(f2cross)
f2result <- scanMQM(f2cross)                    	# Do a MQM scan of the genome
plotMQMone(f2result)					# Use our fancy plotting routine
}

\keyword{ hplot }
