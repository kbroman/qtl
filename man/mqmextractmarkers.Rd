\name{mqmextractmarkers}
\alias{mqmextractmarkers}

\title{MQM marker extraction}

\description{
  Extract the real markers from a cross object that includes pseudo markers
}

\usage{
mqmextractmarkers(mqmresult)
}

\arguments{
  \item{mqmresult}{result from \code{mqmscan}, including pseudo markers}
}

\value{
  Returns a scanone object with the pseudo markers removed
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{broman@wisc.edu} % \mqmauthors
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
data(multitrait)

multitrait <- fill.geno(multitrait)

result <- mqmscan(multitrait)
newresult <- mqmextractmarkers(result)
}

\keyword{utilities}
