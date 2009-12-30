\name{mqmcofactors}
\alias{mqmcofactors}

\title{Create cofactors to be used with mqmscan }

\description{
  Return a cofactorslist formatted to be used with the mqmscan function
}

\usage{
mqmcofactors(cross, cofactors, sexfactors, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
  \item{cofactors}{ List of integer describing which marker should be treated as a cofactor }
  \item{sexfactors}{ List of integer describing which marker should be treated as a sexfactor }
  \item{verbose}{If TRUE, print tracing information.}
}

\value{
An formatted list of cofactors to be used with the \code{\link{mqmscan}} algorithm
}

\author{ 
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
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
qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
a <- mqmcofactors(cross,c(20,30,80,90),c(186,187))	# Set Cofactors on marker 20,30,80 & 90 and sexfactors on 186,187
result <- mqmscan(cross,a)				# Do a backward model selection
}

\keyword{ models }
