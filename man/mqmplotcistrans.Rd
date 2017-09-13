\name{mqmplot.cistrans}
\alias{mqmplot.cistrans}

\title{cis-trans plot}

\description{
Plot results for a genomescan using a multiple-QTL model. With genetic location
for the traits it is possible to show cis- and trans- locations, and detect
trans-bands
}

\usage{
mqmplot.cistrans(result, cross, threshold=5, onlyPEAK=TRUE,
                 highPEAK=FALSE, cisarea=10, pch=22, cex=0.5,
                 verbose=FALSE, \dots)
}

\arguments{
  \item{result}{An object of class \code{mqmmulti}. See \code{\link{mqmscanall}} for details.}
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{threshold}{ Threshold value in LOD, Markers that have a
        LOD score above this threshold are plotted as small squares
        (see \code{pch} parameter). The markers with LODscores below
        this threshold are not visible
  }
  \item{onlyPEAK}{ Plot only the peak markers ? (TRUE/FALSE)
        (Peak markers are markers that have a QTL likelihood above
        \code{threshold} and higher than other markers in the same region)
  }
  \item{highPEAK}{ Highlight peak markers ? (TRUE/FALSE).
        When using this option peak markers (the marker with the
        highest LOD score in a region above the threshold gets
        an 25\% increase in size and is displayed in red)
  }
  \item{cisarea}{ Adjust the two green lines around the line y=x }
  \item{pch}{ What kind of character is used in plotting of the figure (Default: 22, small square) }
  \item{cex}{ Size of the points plotted (default to 0.5 half of the original size)}
  \item{verbose}{If TRUE, give verbose output}
  \item{\dots}{ Extra parameters will be passed to points }
}

%\details{
%
%}

\value{
        Plotting routine, so no return
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
data(multitrait)
\dontshow{multitrait <- subset(multitrait, chr=1,
                               ind=!apply(multitrait$pheno,
                                          1, function(a) any(is.na(a))))}
data(locations)
multiloc <- addloctocross(multitrait,locations)
multiloc <- calc.genoprob(multiloc)
results <- scanall(multiloc, method="hk")
mqmplot.cistrans(results, multiloc, 5, FALSE, TRUE)
}

\keyword{ hplot }
