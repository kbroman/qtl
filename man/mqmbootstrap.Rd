\name{mqmpermutation}
\alias{mqmpermutation}

\title{Estimate QTL LOD score significance using permutations or simulations }

\description{
   Different randomization approaches to obtain estimates of QTL significance, currently two methods are implemented:
   \itemize{
   \item Random redistribution of traits (method=0)
   \item Random redistribution of simulated trait values (method=1)
   }
   Calculations can be parallelized using the SNOW package.
}

\usage{
	mqmpermutation(cross, scanfunction=scanone, pheno.col=1,
	          multicore=TRUE, n.perm=10, b.size=10, file="MQM_output.txt",
   	          n.cluster=1, method=c("permutation","simulation"), plot=FALSE, verbose=FALSE, \dots)
}

\arguments{
	\item{cross}{ 
An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
	\item{scanfunction}{ Function to use when mappingQTL's (either scanone,cim or mqm)}
	\item{pheno.col}{ 
Column number in the phenotype matrix which should be used as the phenotype. This can be a vector of integers. % \phenocol
  }
	\item{multiC}{ Use multicore (if available)}  
	\item{n.run}{ Number of bootstrapping runs (DEFAULT=10, should be 1000, or higher,
  for publications) }
	\item{b.size}{ Batch size. The entire set is split in jobs. Each job contains b.size number of traits per job }
        \item{file}{ Name of the intermediate output file used }
	\item{n.cluster}{ Number of child processes to split the job into }
	\item{method}{ What kind bootstrapping should occur: permutation or simulation } 
	\item{plot}{If TRUE, make a plot}
	\item{verbose}{If TRUE, print tracing information}
	\item{\dots}{Parameters passed through to the
	  \code{\link{scanone}}, \code{\link{cim}} or
	    \code{\link{mqm}} functions} 
}

\details{
  Analysis of \code{\link{scanone}}, \code{\link{cim}} or
	    \code{\link{mqm}} to scan for QTL in shuffled/randomized data. It is recommended to also install the \code{snow} library.
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
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmsetcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmpermutation}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
	}
}

\examples{
data(multitrait)				#Use the multitrait dataset
\dontshow{multitrait <- subset(multitrait, chr=1:2, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
multitrait$pheno <- log(multitrait$pheno)}
multitrait <- calc.genoprob(multitrait)
result <- mqmpermutation(multitrait,pheno.col=7, n.run=2, b.size=2)}
multitrait <- fill.geno(multitrait)
cof <- mqmsetcofactors(multitrait,3)		#Set cofactors at each 3th marker
\dontrun{f2result <- mqmpermutation(multitrait,cofactors=cof,pheno.col=7, b.size=5,verbose=FALSE)}
\dontshow{f2result <- mqmpermutation(multitrait,cofactors=cof,pheno.col=7, n.run=2,b.size=2,verbose=FALSE)}
#EXTRA
f2perm <- mqmprocesspermutation(f2result)			#Create a permutation object
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
