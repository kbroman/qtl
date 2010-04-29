\name{simulateMissingData}
\alias{simulateMissingData}

\title{Simulates missing genotype data}

\description{
  Simulate missing genotype data by randomly dropping out genotype marker annotations
}

\usage{
  simulateMissingData(cross, percentage = 5)
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
multimissing5 <- simulateMissingData(multitrait,perc=5)
perc <- (sum(nmissing(multimissing5))/sum(ntyped(multimissing5)))
}

\keyword{utilities}
