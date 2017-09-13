\name{mqmplot.singletrait}
\alias{mqmplot.singletrait}

\title{Plot LOD curves of a multiple-QTL model}

\description{
  Plot the LOD curve for a genome scan for a single trait, with a multiple-QTL model (the
  output of \code{\link{mqmscan}}).
}

\usage{
mqmplot.singletrait(result, extended = 0 ,\dots)
}

\arguments{
  \item{result}{ \code{\link{mqmscan}} result.}
  \item{extended}{ Extended plotting of the information content }
  \item{\dots}{ Extra arguments passed to \code{\link{plot.scanone}} }
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
#Simulated F2 Population
f2qtl <- c(3,15,1,0)                                    # QTL at chromosome 3
data(map10)                                             # Mouse genetic map
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
f2cross <- mqmaugment(f2cross)
f2result <- mqmscan(f2cross)                        # Do a MQM scan of the genome
mqmplot.singletrait(f2result)					# Use our fancy plotting routine
}

\keyword{ hplot }
