\name{mqmscanfdr}
\alias{mqmscanfdr}

\title{ Estimate FDR for multiple trait QTL analysis }

\description{
  Estimate the false discovery rate (FDR) for multiple trait analysis
}

\usage{
mqmscanfdr(cross, scanfunction=mqmscanall, 
  thresholds=c(1,2,3,4,5,7,10,15,20), n.perm=10, 
  verbose=FALSE, \dots
  )
}

\arguments{
	\item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
	\item{scanfunction}{ QTL mapping function, Note: Must use scanall or mqmscanall. Otherwise this will 
  not produce usefull results. Reason: We need a function that maps all traits ecause of the correlation
  structure which is not changed (between traits) during permutation (Valis options: scanall or mqmscanall)
  }
	\item{thresholds}{ False discovery rate (FDR) is calculated for peaks 
    above these LOD thresholds (DEFAULT=Range from 1 to 20, using 10 thresholds)
    Parameter is a list of LOD scores at which FDR is calculated.
  }
	\item{n.perm}{ Number of permutations (DEFAULT=10 for quick analysis, however
    for publications use 1000, or higher)
  }
	\item{verbose}{ verbose output }
	\item{\dots}{ Parameters passed to the mapping function
  } 
}

\details{
  This function wraps the analysis of \code{\link{scanone}}, \code{\link{cim}}
  and \code{\link{mqmscan}} to scan for QTL in shuffled/randomized data. It is 
  recommended to also install the \code{snow} library for parallelization of
  calculations. The \code{snow} library allows
  calculations to run on multiple cores or even scale it up to an entire cluster, 
  thus speeding up calculation by the number of computers used.
}

\value{
	Returns a data.frame with 3 columns: FalsePositives, FalseNegatives and False Discovery Rates. 
  In the rows the userspecified thresholds are with scores for the 3 columns.
}

\author{ 
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
  }

\seealso{
	\itemize{
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
data(multitrait)

multitrait <- fill.geno(multitrait) # impute missing genotype data
\dontrun{  result <- mqmscanfdr(multitrait, threshold=10.0, n.perm=1000) # Calculate the thresholds
}\dontshow{
  multitrait <- subset(multitrait, chr=1, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
  multitrait$pheno <- multitrait$pheno[,1:3]
  result <- mqmscanfdr(multitrait, threshold=10.0, n.perm=2)
}
}

\references{
\itemize{
% \input{"inst/docs/Sources/MQM/mqm/significance_references.txt"}
    \item Bruno M. Tesson, Ritsert C. Jansen (2009) Chapter 3.7. Determining the significance threshold \emph{eQTL Analysis in Mice and Rats} \bold{1}, 20--25
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963--971.
% -----^^ inst/docs/Sources/MQM/mqm/significance_references.txt ^^-----
% \input{"inst/docs/Sources/MQM/mqm/parallelisation_references.txt"}
		\item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
		\item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.2-1. 	
% -----^^ inst/docs/Sources/MQM/mqm/parallelisation_references.txt ^^-----
}
}

\keyword{ models }
