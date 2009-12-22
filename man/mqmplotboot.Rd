\name{mqmplot_permutations}
\alias{mqmplot_permutations}
\alias{mqmplot_boot}

\title{Plot results from a bootstrap scan}

\description{
	Plotting routine to display the results from a bootstrap scan 
}

\usage{
mqmplot_permutations(result, \dots)
}

\arguments{
  \item{result}{ Resulting \code{mqmmulti} object created when performing 
	permutation analysis. Only permutations done by the MQM algorithm \code{\link{mqm}} 
	can use this plotting routine.
  }
  \item{\dots}{ Extra argument passed to the polyplot method of R. Brouwer }    
}

%\details{
%  ~~ If necessary, more details than the description above ~~
%}

\value{ No value returned (plotting routine) }

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; Rutger Brouwer \email{R.W.W.Brouwer@rug.nl} }

\seealso{ 
\itemize{
	\item \code{\link{bootstrap}} - Function called to estimate confidence interfalls for scanone, cim and mqm
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
mqmplot_boot(f2result)					# Use the nice bootstrap plotting methodes
}
\keyword{ hplot }
