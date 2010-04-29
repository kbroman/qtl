\name{mqmsetcofactors}
\alias{mqmsetcofactors}

\title{Create cofactors to be used with mqm} 

\description{
  Return a list of cofactors ready to be used with the mqm function
}

\usage{
mqmsetcofactors(cross, each = NULL, cofactors=NULL, sexfactors=NULL, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.} 
  \item{each}{ Every 'each' marker will be used as a cofactor, when each is used the \code{cofactors} and \code{sexfactors} parameter is ignored }
  \item{cofactors}{ List of integer describing which markers should be treated as additive cofactors, when \code{each} is set, this parameter is ignored }
  \item{sexfactors}{ List of integer describing which markers should be treated as dominant cofactors (sexfactors), when \code{each} is set, this parameter is ignored }
  \item{verbose}{If TRUE, print tracing information.}
}

\value{
	An formatted list of cofactors to be used with the $cofactors$ parameter of the \code{\link{MQM}} algorithm.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{ 
  \itemize{
    % \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmsetcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
	data(hyper)							        # Hyper dataset
	\dontshow{hyper <- hyper[1:5]}
	hyperfilled <- fill.geno(hyper)
  # Automatic cofactors every third marker
	cofactors <- mqmsetcofactors(hyperfilled,3)
	result <- mqmscan(hyperfilled,cofactors)	# Backward model selection
	mqmgetmodel(result)
  #Manual cofactors at markers 3,6,9,12,40 and 60
  cofactors <- mqmsetcofactors(hyperfilled,cofactors=c(3,6,9,12,40,60))	
	result <- mqmscan(hyperfilled,cofactors)	# Backward model selection
	mqmgetmodel(result)
}

\keyword{ models }
