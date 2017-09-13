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
  The function returns the multiple QTL model created, which consists of the cofactors selected during the modeling phase of the algorithm.
  This model was used when scanning for additional QTL in the mqmscan function. The format of the model is compatible with the
  \code{\link{makeqtl}} function. For more information about the format of the model see the \code{\link{makeqtl}} page.
  When no cofactor was selected in the modeling phase no model was created, then this function will return a NULL value.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{
    \itemize{
     \item \code{\link{mqmsetcofactors}} - Setting multiple cofactors for backward elimination
    \item \code{\link{makeqtl}} - Make a qtl object
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
    data(hyper)
        \dontshow{hyper <- subset(hyper, chr=1:4)}
    hyperfilled <- fill.geno(hyper)
    cofactors <- mqmsetcofactors(hyperfilled,4)
    result <- mqmscan(hyperfilled,cofactors)
    mqmgetmodel(result)
    plot(mqmgetmodel(result))
}

\keyword{ utilities }
