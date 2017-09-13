\name{mqmplot.cofactors}
\alias{mqmplot.cofactors}

\title{Plot cofactors on the genetic map}

\description{
  Plots cofactors as created by \code{\link{mqmsetcofactors}}
  or \code{\link{mqmautocofactors}} on the genetic map.
}

\usage{
mqmplot.cofactors(cross,cofactors, \dots)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{cofactors}{
List of cofactors to be analysed in the QTL model. To set cofactors use \code{\link{mqmautocofactors}} or \code{mqmsetcofactors}}. % \mqmcofactors
  \item{\dots}{Passed to \code{\link{plot.qtl}}}
}

\author{
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
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
data(multitrait)
cof1 <- mqmsetcofactors(multitrait,20)
cof2 <- mqmsetcofactors(multitrait,10)
op <- par(mfrow=c(2,1))
mqmplot.cofactors(multitrait,cof1,col="blue")
mqmplot.cofactors(multitrait,cof2,col="blue")
op <- par(mfrow=c(1,1))
}

\keyword{ hplot }
