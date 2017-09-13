\name{mqmscanall}
\alias{mqmscanall}

\title{Parallelized MQM on multiple phenotypes in a cross object}

\description{
  Parallelized QTL analysis using MQM on multiple phenotypes in a cross object (uses SNOW)
}

\usage{
    mqmscanall(cross, multicore=TRUE, n.clusters = 1,batchsize=10,cofactors=NULL, \dots)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{multicore}{ Use multiple cores (only if the package SNOW is available, otherwise this setting will be ignored)
  }
  \item{n.clusters}{ Number of parallel processes to spawn, recommended is setting this lower than the number of cores in the computer
  }
  \item{batchsize}{ Batch size. The entire set is split in
    jobs to reduce memory load per core. Each job contains batchsize number of traits per job.
  }
  \item{cofactors}{cofactors, only used when scanfunction is mqmscan.
List of cofactors to be analysed in the QTL model. To set cofactors use \code{\link{mqmautocofactors}} or \code{mqmsetcofactors}}. % \mqmcofactors
  \item{\dots}{ Parameters passed through to the \code{\link{mqmscan}} function used in
    scanning for QTLs
  }
}

\details{
  Uses \code{mqmscan} to scan for QTL's for each phenotype in the cross object.
  It is recomended that the package SNOW is installed before using this function
  on large numbers of phenotypes.
}

\value{
  Returns a MQMmulti object. This object is a list of \code{\link{scanone}}
  objects that can be plotted using \code{plot.scanone(result[[trait]])} or using \code{mqmplot.multitrait(result)}
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
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
#Doing a multitrait analysis
data(multitrait)
\dontshow{multitrait <- subset(multitrait, chr=1:2, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
multitrait$pheno <- multitrait$pheno[,1:3]}
multitrait <- calc.genoprob(multitrait)
cof <- mqmsetcofactors(multitrait,3)
multitrait <- fill.geno(multitrait)
result <- mqmscanall(multitrait,cofactors=cof,batchsize=5)
mqmplot.multitrait(result,"lines")
}

\references{
    \itemize{
% \input{"inst/doc/Sources/MQM/mqm/parallelisation_references.txt"}
        \item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
        \item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.2-1.
% -----^^ inst/doc/Sources/MQM/mqm/parallelisation_references.txt ^^-----
    }
}

\keyword{ models }
