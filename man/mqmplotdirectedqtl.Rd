\name{mqmplot_directedqtl}
\alias{mqmplot_directedqtl}

\title{Plotting routine to display the results from a QTL scan  }

\description{
  Plotting routine to display the results from a mqmscan
}

\usage{
mqmplot_directedqtl(cross, scanone, pheno.col=1, draw = TRUE)
}

\arguments{
  \item{cross}{
    An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
  \item{result}{ Result object from scanone or mqmscan }
  \item{pheno.col}{ Which phenotype column was originaly scanned for 
                    QTLs, the pheno.col specified here will be used 
                    to extract QTL effects }
  \item{draw}{ Should the resulting profile be plotted ? setting draw to FALSE 
               will just return the scanone object with LOD scores }
}

\author{ 
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{ 
  \itemize{
    % \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmsetcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
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
f2resultd <- mqmplot_directedqtl(f2cross,f2result)					# Use our fancy plotting routine
}

\keyword{ hplot }
