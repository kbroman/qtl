\name{mqmplot.directedqtl}
\alias{mqmplot.directedqtl}

\title{Plot LOD*Effect curves of a multiple-QTL model }

\description{
  Plot the LOD*Effect curve for a genome scan with a multiple-QTL model (the
  output of \code{\link{mqmscan}}).
}

\usage{
mqmplot.directedqtl(cross, mqmresult, pheno.col=1, draw = TRUE)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{mqmresult}{ Results from mqmscan of type \code{scanone} }
  \item{pheno.col}{ From which phenotype in the crossobject are the result calculated }
  \item{draw}{ If TRUE, draw the figure. }
}

\value{
  Returns a scanone object, with added the effectsign calculated
  internally by the function effect.scan. For more info on the
  scanone object see: \code{\link{scanone}}
}

\author{
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{
  \itemize{
    % \input{"inst/doc/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{http://rqtl.org/tutorials/MQM-tour.pdf}
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
f2cross <- fill.geno(f2cross)                           # Fill in missing genotypes
f2result <- mqmscan(f2cross)                        # Do a MQM scan of the genome
mqmplot.directedqtl(f2cross,f2result)
}

\keyword{ hplot }
