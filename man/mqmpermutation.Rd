\name{mqmpermutation}
\alias{mqmpermutation}

\title{Estimate QTL LOD score significance using permutations or simulations }

\description{
   Two randomization approaches to obtain estimates of QTL significance:
   \itemize{
   \item Random redistribution of traits (method='permutation')
   \item Random redistribution of simulated trait values (method='simulation')
   }
   Calculations can be parallelized using the SNOW package.
}

\usage{
	mqmpermutation(cross, scanfunction=scanone, pheno.col=1,
	          multicore=TRUE, n.perm=10, batchsize=10, file="MQM_output.txt",
   	          n.cluster=1, method=c("permutation","simulation"),cofactors=NULL, plot=FALSE, verbose=FALSE, \dots)
}

\arguments{
	\item{cross}{ 
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
	\item{scanfunction}{ Function to use when mappingQTL's (either scanone,cim or mqm)}
	\item{pheno.col}{ 
Column number in the phenotype matrix which should be used as the phenotype. This can be a vector of integers. % \phenocol
  }
	\item{multicore}{ Use multicore (if available)}  
	\item{n.perm}{ Number of permutations to perform (DEFAULT=10, should be 1000, or higher,
  for publications) }
	\item{batchsize}{ Batch size. The entire set is split in jobs. Each job contains b.size number of traits per job }
        \item{file}{ Name of the intermediate output file used }
	\item{n.cluster}{ Number of child processes to split the job into }
	\item{method}{ What kind permutation should occur: permutation or simulation } 
  \item{cofactors}{cofactors, only used when scanfunction is mqm.
List of cofactors to be analysed in the QTL model. To set cofactors use \code{\link{mqmautocofactors}} or \code{mqmsetcofactors}}. % \mqmcofactors
	\item{plot}{If TRUE, make a plot}
	\item{verbose}{If TRUE, print tracing information}
	\item{\dots}{Parameters passed through to the
	  \code{\link{scanone}}, \code{\link{cim}} or
	    \code{\link{mqmscan}} functions} 
}

\details{
  Analysis of \code{\link{scanone}}, \code{\link{cim}} or
	    \code{\link{mqmscan}} to scan for QTL in shuffled/randomized data. It is recommended to also install the \code{snow} library.
  The \code{snow} library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}

\value{
	Returns a mqmmulti object. this object is a list of scanone objects that can be plotted using \code{plot.scanone(result[[trait]])}
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
data(multitrait)				#Use the multitrait dataset
\dontshow{multitrait <- subset(multitrait, chr=1:2, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
multitrait$pheno <- log(multitrait$pheno)}
multitrait <- calc.genoprob(multitrait)
result <- mqmpermutation(multitrait,pheno.col=7, n.perm=2, batchsize=2)

\dontrun{cof <- mqmautocofactors(multitrait,50)		#Set 50 cofactors
}\dontshow{cof <- mqmautocofactors(multitrait,5)}
multitrait <- fill.geno(multitrait)
result <- mqmpermutation(multitrait,scanfunction=mqmscan,cofactors=cof,pheno.col=7, n.perm=2,batchsize=2,verbose=FALSE)
#EXTRA
f2perm <- mqmprocesspermutation(result)			#Create a permutation object
summary(f2perm)						#Get Significant LOD thresholds
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
