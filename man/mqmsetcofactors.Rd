\name{mqmsetcofactors}
\alias{mqmsetcofactors}

\title{Set cofactors at fixed distances, to be used with MQM} 

\description{
  Sets cofactors, at fixed distance. Together
  with \code{mqmscan} cofactors are selected through backward elimination.
}

\usage{
mqmsetcofactors(cross, each = NULL, cofactors=NULL, sexfactors=NULL, verbose=FALSE)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
   } 
  \item{each}{ Use every marker in list as a cofactor }
  \item{cofactors}{ 
List of cofactors to be analysed as cofactors in the QTL model. To set cofactors use \code{\link{mqmautocofactors}} or \code{mqmsetcofactors}}. % \mqmcofactors
  \item{sexfactors}{ List of integer describing which marker should be treated as a sexfactor }
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
    \item MQM-tutorial - The online MQM tutorial
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Find the most informative cofactors through backward elimination
    \item \code{\link{mqmsetcofactors}} - Set cofactors
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
	data(hyper)							        # Hyper dataset
	\dontshow{hyper <- hyper[1:5]}
	hyperfilled <- fill.geno(hyper)
	cofactors <- mqmsetcofactors(hyperfilled,3)	# Set Cofactors every third marker
	result <- mqmscan(hyperfilled,cofactors)	# Backward model selection
	mqmgetmodel(result)
}

\keyword{ models }
