\name{mqmgetmodel}
\alias{mqmgetmodel}
\alias{mqmgetmodel}

\title{ Retrieve the QTL model used in mapping from the results of an MQM scan }

\description{
  Retrieves the QTL model used for scanning from the output of an MQM scan. The
  model only contains the selected cofactors significant at the specified
  cofactor.significance from the results of an mqm scan
}

\usage{
	mqmgetmodel(scanresult)
}

\arguments{
  \item{scanresult}{
An object returned by \code{mqmscan}, including cofactors and QTL model. % \mqmscanobject
  } 
}

\value{
  The function returns the QTL model created (cofactors selected in modeling
  found to be significant at cofactor.significance) and used during the mqmscan
  See for the format of the model: \code{\link{makeqtl}} returns NULL when no
  cofactors reach the significance threshold during QTL modeling.
}  

\author{ 
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{
	\itemize{
     \item \code{\link{mqmsetcofactors}} - Setting multiple cofactors for backward elimination
    \item \code{\link{makeqtl}} - Make a qtl object
   % \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{http://www.rqtl.org/tutorials/MQM-tour.pdf} 
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Set cofactors using marker density
    \item \code{\link{mqmsetcofactors}} - Set cofactors at fixed locations
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
	}
}

\examples{
	data(hyper)
        \dontshow{hyper <- subset(hyper, chr=1:4)}
	hyperfilled <- fill.geno(hyper)
	cofactors <- mqmsetcofactors(hyperfilled,4)
	result <- mqmscan(hyperfilled,cofactors)
	mqmgetmodel(result)
	plot(mqmgetmodel(result))
}

\keyword{ utilities }
