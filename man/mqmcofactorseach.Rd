\name{mqmsetcofactors}
\alias{mqmsetcofactors}

\title{Create cofactors to be used with mqm} 

\description{
  Return a list of cofactors ready to be used with the mqm function
}

\usage{
mqmsetcofactors(cross, each = 3, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.} 
  \item{each}{ Every 'each' marker will be used as a cofactor }
  \item{verbose}{If TRUE, print tracing information.}
}

\value{
	An formatted list of cofactors to be used with the $cofactors$ parameter of the \code{\link{mqm}} algorithm.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{ 
  \itemize{
    % \input{"inst/docs/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
data(map10)							                  # Mouse genome
qtl <- c(3,15,3,7)						            # QTL at chromosome 3
\dontshow{map10 <- map10[1:5]}
cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
a <- mqmsetcofactors(cross,3)					    # Set Cofactors every third marker
result <- mqmscan(cross,a)					      # Backward model selection
}

\keyword{ models }
