\name{mqmplot.cofactors}
\alias{mqmplot.cofactors}

\title{Plotting routine to display cofactors on the genetic map}

\description{
  Plotting routine to display cofactors created by \code{\link{mqmautocofactors}} 
  or \code{\link{mqmsetcofactors}} on the genetic map. 
}

\usage{
mqmplot.cofactors(cross,cofactors, \dots)
}

\arguments{
  \item{cross}{ 
    An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
  \item{cofactors}{
    List of cofactors to be analysed as cofactors in the 
    QTL model. See \code{\link{mqmsetcofactors}} on how-to set cofactors 
    for backward elimination. (0:no cofactor, 1:cofactor, 2:sexfactor) 
  }
  \item{\dots}{Passed to \code{\link{plot.qtl}}}
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
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
data(multitrait)
autocof <- mqmautocofactors(multitrait,20)
eachcof <- mqmsetcofactors(multitrait,10)
op <- par(mfrow=c(2,1))
mqmplot.cofactors(multitrait,autocof, col="blue")
mqmplot.cofactors(multitrait,eachcof, col="blue")
op <- par(mfrow=c(1,1))
}

\keyword{ hplot }
