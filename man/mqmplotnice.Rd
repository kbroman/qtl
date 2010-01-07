\name{mqmplot_nice}
\alias{mqmplot_nice}

\title{Plot the results from a mqm scan on multiple phenotypes  }

\description{
	Plotting routine to display the results from a mqm scan on multiple phenotypes
}
\usage{
mqmplot_nice(result, \dots)
}

\arguments{
  \item{result}{ Result object from mqm scan \code{\link{mqm}}. 
    Make sure that multiple phenotypes have been scanned. 
  }
  \item{\dots}{ Extra argument passed to the polyplot methode of 
    R. Brouwer. Use \code{legendloc=1} to plot legend outside of
    plotting area 
  }   
}

\author{ 
Danny Arends \email{danny.arends@gmail.com} ; Rutger Brouwer % \dannyrutgerauthors
}

\seealso{ 
  \itemize{
    % \input{"inst/docs/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/mqm/standard_seealso.txt ^^-----
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
mqmplot_nice(f2result)					# Use the nice plotting methodes
}
\keyword{ hplot }
