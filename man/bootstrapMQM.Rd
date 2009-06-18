\name{bootstrap}
\alias{bootstrap}
\alias{bootstrapcim}
\alias{bootstrapmqm}
\alias{permutate}

\title{randomization approaches to obtain estimates of QTL LOD score significance }

\description{
   Using different randomization approaches to obtain estimates of QTL significance, currently two methods are implemented:
   \itemize{
   \item Random redistribution of traits (method=0)
   \item Random redistribution of simulated trait values (method=1)
   }
}

\usage{
	bootstrap(cross = NULL, Funktie=scanone, pheno.col=1,
	          multiC=TRUE, n.run=10, b_size=10, file="MQM_output.txt",
   	          n.clusters=2, bootmethod=0, plot=FALSE, verbose=FALSE, \dots)
	bootstrapcim(\dots)
	bootstrapmqm(\dots)
}

\arguments{
	\item{cross}{ Cross object, An object of class \code{cross} from R/QTL with multiple phenotypes suitable for scanone and scanMQM from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
	\item{Funktie}{ Function to use when mappingQTL's (either scanone,cim or mqm)}
	\item{pheno.col}{ Phenotype under investigation (from the cross object) }
	\item{multiC}{ Use multicore (if available)}  
	\item{n.run}{ Number of bootstrapping runs (DEFAULT=10, should be 1000 to 5000) }
	\item{b_size}{ Batch size. The entire set is split in jobs. Each job contains b_size # of traits per job }
        \item{file}{File name (?)}
	\item{n.clusters}{ Number of child processes to split the job into }
	\item{bootmethod}{ What kind of bootstrapping method should be used ( 0= Random Shuffle (DEFAULT)/ 1= Parametric ) } 
	\item{plot}{If TRUE, make a plot}
	\item{verbose}{If TRUE, print tracing information}
	\item{\dots}{Parameters passed through to the
	  \code{\link[qtl]{scanone}}, \code{\link[qtl]{cim}} or
	    \code{\link[qtl]{mqm}} functions} 
}

\details{
  Analysis of \code{\link[qtl]{scanone}}, \code{\link[qtl]{cim}} or
	    \code{\link[qtl]{mqm}} to scan for QTL in shuffled/randomized data. It is recommended to also install the \code{snow} library.
  The \code{snow} library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}

\value{
	Returns a MQMmulti object. this object is a list of scanone objects that can be plotted using \code{plot.scanone(result[[trait]])}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} ; Karl W. Broman \email{kbroman@biostat.wisc.edu} }

\note{ 
	If you find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{
	\itemize{
	\item \code{\link[qtl]{MQMCofactors}} - Learn how to set Cofactors for the MQM algorithm
	\item \code{\link[qtl]{scanall}} - Function called to analyse multiple traits in the cross object
	\item \code{\link[qtl]{scanone}} - Function called to do single trait analysis
	\item \code{\link[qtl]{cim}} - Function called to do single trait analysis
	\item \code{\link[qtl]{mqm}} - Function called to do single trait analysis
	}
}

\examples{
data(multitrait)				#Use the multitrait dataset
result <- bootstrap(multitrait,pheno.col=7, b_size=5)
result <- bootstrapcim(multitrait,pheno.col=7, b_size=5)
cof <- MQMCofactorsEach(multitrait,3)		#Set cofactors at each 3th marker
f2result <- bootstrapmqm(multitrait,cofactors=cof,pheno.col=7, b_size=5,verbose=FALSE)	# Bootstrap 200 runs in batches of 20

#EXTRA
f2perm <- MQMpermObject(f2result)			#Create a permutation object
summary(f2perm)						#Get Significant LOD thresholds

}

\references{
	\itemize{
		\item Jansen RC, (2007) Quantitative trait loci in inbred lines. Chapter 18 of \emph{Handbook of Stat. Genetics} 3rd edition. John Wiley & Sons, Ltd.
		\item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.2-1. 	
		\item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388--391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 1447--1455.
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963--971.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205--211.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 1--38.
	}
}

\keyword{ models }
