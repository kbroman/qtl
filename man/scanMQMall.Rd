\name{scanall}
\alias{scanall}
\alias{mqmall}
\alias{cimall}
\title{ 
	scanall - Multitrait analysis of a cross object 
}

\description{
	Uses scanone, cim, or mqm to scan for QTL's for each phenotype in the cross object.
}
\usage{
	scanall(cross = NULL, Funktie=scanone, n.clusters = 2,multiC=TRUE, b_size=10, plot=TRUE, verbose=TRUE,\dots)
	cimall(cross = NULL,\dots)
	mqmall(cross = NULL,\dots)
}

\arguments{
  \item{cross}{ Cross object, an object of class \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
  \item{Funktie}{ Function to use when mappingQTL's (either scanone,cim, mqm or your own mapping function). Easier call the interface implementing the function under interest} 
  \item{multiC}{ Use multicore (if available)}
  \item{n.clusters}{ Number of child processes to split the job into }
  \item{b_size}{ Batch size. The entire set is split in jobs. Each job contains b_size # of traits per job }
  \item{\dots}{Parameters passed through to the scanone, cim or mqm functions}
}

\details{
  Uses the function specified by the Funktie parameter to scan for QTL's for each phenotype in the cross object. It is recomended to also install the snow library.
  The snow library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}

\value{
  Returns a MQMmulti object. This object is a list of scanone objects that can be plotted using \code{plot.scanone(result[[trait]])} or using \code{plot.all(result)}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu}  }

\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{ 
	\itemize{
	\item \code{\link{multitrait}} - Dataset used in the example analysis
	\item \code{\link{fill.geno}} - Augmentationroutine to estimate missing data	
	\item \code{\link{MQMaugment}} - Augmentationroutine to estimate missing data
	\item \code{\link{scanone}} - Function called to do single trait analysis
	\item \code{\link{cim}} - Function called to do single trait analysis
	\item \code{\link{mqm}} - Function called to do single trait analysis
	}
}

\examples{
	#Doing a multitrait analysis
	data(multitrait)
	result <- scanall(multitrait, b_size=5)
	result <- cimall(multitrait, b_size=5)
	cof <- MQMCofactorsEach(multitrait,3)
	result <- mqmall(multitrait,cofactors=cof,doLOG=1,b_size=5)

	#FUNNY, run at own risk: 3D QTL screensaver
	phi<-0
	while(T){
		plot.all(result,"D",phi);phi<- phi+10;Sys.sleep(0.01)
	}
}
\references{
	\itemize{
		\item Jansen RC, (2007), Chapter 18 - Quantitative trait loci in inbred lines. \emph{Handbook of Stat. Genetics} 3th edition,(c) 2007 John Wiley & Sons, Ltd.
		\item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.2-1. 	
		\item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388-391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 1447-1455.
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963--971.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205--211.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 1--38.
	}
}
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
