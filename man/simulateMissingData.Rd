\name{simulatemissingdata}
\alias{simulatemissingdata}

\title{Simulates missing genotype data}

\description{
  Simulate missing genotype data by removing some genotype data from the cross object
}

\usage{
  simulatemissingdata(cross, percentage = 5)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{percentage}{ How much of the genotype data do we need to randomly drop? }
}

\author{
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\value{
  An object of class \code{cross} with \code{percentage} % of (randomly) missing genotype data. See \code{\link{read.cross}} for details
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
multimissing5 <- simulatemissingdata(multitrait,perc=5)
perc <- (sum(nmissing(multimissing5))/sum(ntyped(multimissing5)))
}

\keyword{utilities}
