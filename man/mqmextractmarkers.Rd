\name{mqmextractmarkers}
\alias{mqmextractmarkers}

\title{MQM marker extraction}

\description{
  Extract only the real markers from an cross object
}

\usage{
mqmextractmarkers(mqmresult)
}

\arguments{
  \item{mqmresult}{result from \code{mqmscan} to scale down to only the results at the real markers}
}

\value{
  Returns a scanone object with the pseudomarkers removed, containing only the real markers.
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
data(multitrait)

multitrait <- fill.geno(multitrait)

result <- mqmscan(multitrait)
newresult <- mqmextractmarkers(result)
}

\keyword{utilities}
