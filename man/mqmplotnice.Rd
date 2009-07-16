\name{mqmplotnice}
\alias{mqmplotnice}

\title{Plot the results from a bootstrap or scanall scan  }

\description{
	Plotting routine to display the results from a bootstrap or scanall scan
}
\usage{
mqmplotnice(result, \dots)
}

\arguments{
  \item{result}{ Result object from scanall }
  \item{\dots}{ Extra argument passed to the polyplot methode of R. Brouwer. Use \code{legendloc=1} to plot legend outside of plotting area }   
}

%\details{
%  ~~ If necessary, more details than the description above ~~
%}

\value{None.}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; Rutger Brouwer \email{R.W.W.Brouwer@rug.nl} }

\seealso{ 
\itemize{
	\item \code{\link[qtl]{scanall}} - Function called to analyse multiple traits
	\item \code{\link[qtl]{bootstrap}} - Function called to estimate confidence interfalls for MQM
}
}

\examples{
#Simulated F2 Population
f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
data(map10)                                             # Mouse genome
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
f2cross <- fill.geno(f2cross)
\dontrun{f2result <- mqmbootstrap(f2cross)                       # Bootstrap MQM to obtain confidence intervals
}\dontshow{f2result <- mqmbootstrap(f2cross, n.run=2, b.size=2)}
mqmplotnice(f2result)					# Use the nice plotting methodes
}
\keyword{ hplot }
