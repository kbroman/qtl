\name{mqmsetcofactors}
\alias{mqmsetcofactors}

\title{Set cofactors at fixed intervals, to be used with MQM}

\description{
  Set cofactors, at fixed marker intervals. Together
  with \code{mqmscan} cofactors are selected through backward elimination.
}

\usage{
mqmsetcofactors(cross, each = NULL, cofactors=NULL, sexfactors=NULL, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
  \item{each}{ Every 'each' marker will be used as a cofactor, when each is used the \code{cofactors} and \code{sexfactors} parameter is ignored }
  \item{cofactors}{
List of cofactors to be analysed in the QTL model. To set cofactors use \code{\link{mqmautocofactors}} or \code{mqmsetcofactors}; when \code{each} is set, this parameter is ignored }
  \item{sexfactors}{ list of markers which should be treated as dominant cofactors (sexfactors), when \code{each} is set, this parameter is ignored }
  \item{verbose}{If TRUE, print tracing information.}
}

\value{
    An list of cofactors to be passed into \code{\link{mqmscan}}.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{
  \itemize{
    % \input{"inst/doc/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{https://rqtl.org/tutorials/MQM-tour.pdf}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Set cofactors using marker density
    \item \code{\link{mqmsetcofactors}} - Set cofactors at fixed locations
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/doc/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
    data(hyper)                                 # Hyper dataset
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
