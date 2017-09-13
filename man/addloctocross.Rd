\name{addloctocross}
\alias{addloctocross}

\title{Add phenotype location into a cross object}

\description{
    Add phenotype location(s) into a cross object
  (with eQTL/pQTL studies)
}

\usage{
    addloctocross(cross, locations=NULL, locfile="locations.txt", verbose=FALSE)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{locations}{R variable holding location information}
  \item{locfile}{ load from a file, see the details section for the layout of the file.}
  \item{verbose}{If TRUE, give verbose output}
}

\details{
inputfile layout:
Num	Name	Chr	cM
1	X3.Hydroxypropyl	4	50.0
Num is the number of the phenotype in the cross object
Name is the name of the phenotype (will be checked against the name already in the cross object at position num
Chr Chromosome
cM	postion from start of chromosome in cM
}

\value{
  The input cross object, with the locations added as an aditional component \code{locations}
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{
    \itemize{
    \item \code{\link{mqmplot.cistrans}} - Cis/trans plot
% \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
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
% See the help file for mqmplot.cistrans
\dontrun{
    data(multitrait)
    data(locations)
    multiloc <- addloctocross(multitrait,locations)
    results <- scanall(multiloc)
    mqmplot.cistrans(results, multiloc, 5, FALSE, TRUE)
  }
}

\keyword{ utilities }
