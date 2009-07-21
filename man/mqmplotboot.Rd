\name{mqmplotboot}
\alias{mqmplotboot}

\title{Plot results from a bootstrap scan}

\description{
	Plotting routine to display the results from a bootstrap scan 
}

\usage{
mqmplotboot(result, \dots)
}

\arguments{
  \item{result}{ Result object from scanall }
  \item{\dots}{ Extra argument passed to the polyplot methode of R. Brouwer }    
}

%\details{
%  ~~ If necessary, more details than the description above ~~
%}

\value{None.}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; Rutger Brouwer \email{R.W.W.Brouwer@rug.nl} }

\seealso{ 
\itemize{
	\item \code{\link[qtl]{bootstrap}} - Function called to estimate confidence interfalls for scanone, cim and mqm
}
}

\examples{
#Simulated F2 Population
f2qtl <- c(3,15,3,7)                                    # QTL at chromosome 3
data(map10)                                             # Mouse genome
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
f2cross <- calc.genoprob(f2cross)
\dontrun{f2result <- bootstrap(f2cross)               	        # Bootstrap to obtain confidence intervals
}\dontshow{f2result <- bootstrap(f2cross, n.run=2, b.size=2, method="hk")}
mqmplotboot(f2result)					# Use the nice bootstrap plotting methodes
}
\keyword{ hplot }
