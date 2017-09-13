\name{mqmautocofactors}
\alias{mqmautocofactors}

\title{Automatic setting of cofactors, taking marker density into account}

\description{
  Sets cofactors, taking underlying marker density into account. Together
  with \code{mqmscan} cofactors are selected through backward elimination.
}

\usage{
mqmautocofactors(cross, num=50, distance=5, dominance=FALSE, plot=FALSE, verbose=FALSE)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{num}{ Number of cofactors to set (warns when setting too many cofactors).}
  \item{distance}{ Minimal distance between two cofactors, in cM.}
  \item{dominance}{If TRUE, create a cofactor list that is safe to use
  with the dominance scan mode of MQM. See \code{\link{mqmscan}} for details.}
  \item{plot}{If TRUE, plots a genetic map displaying the selected markers as cofactors.}
  \item{verbose}{If TRUE, give verbose output.}
}

\value{
    A list of cofactors to be used with \code{\link{mqmscan}}.
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
    data(hyper)                     # hyper dataset
    \dontshow{hyper <- hyper[1:5]}
    hyperfilled <- fill.geno(hyper)
    cofactors <- mqmautocofactors(hyperfilled,15)	# Set 15 Cofactors
    result <- mqmscan(hyperfilled,cofactors)	# Backward model selection
    mqmgetmodel(result)
}

\keyword{ models }
