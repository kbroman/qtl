\name{scanMQMall}
\alias{scanMQMall}

\title{ 
	scanMQMall - Multitrait analysis of a cross object 
}

\description{
	Uses MQMaugment and scanMQM to scan for QTL's for each phenotype in the cross object.
}
\usage{
	scanMQMall(cross = NULL, cofactors = NULL, doLOG = 0,REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.clusters = 10,doLOG=0,est.map=0,dominance=0,FF=0)
}

\arguments{
  \item{cross}{ Cross object, an object of class \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{doLOG}{ Indicate if the algorithm can LOG transform traits with variation > 10.000 (DEFAULT 0 = off | 1 = on) }  
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{doLOG}{ Indicate if the algorithm should LOG transform traitdata (DEFAULT 0 = off | 1 = on) }    
  \item{est.map}{ Reestimation of the marker locations on the genetic map before mapping QTLs (no map is returned the old map remains in the cross object }
  \item{dominance}{ Should dominance be taken into account ( 0= Additive  (DEFAULT)/ 1= Dominance) }  
  \item{FF}{ Print an output file at the end of scanning in Frank.Format (file="frank.out") (DEFAULT 0 = no | 1 = yes) }  
  \item{forceRIL}{ Should we force the cross as a RIL, this option removed AB genotypes and replaces them with BB ? ( 0= NO  (DEFAULT)/ 1= YES) }  
  \item{n.clusters}{ Number of child processes to split the job into }
  \item{b_size}{ Batch size. The entire set is split in jobs. Each job contains b_size # of traits per job }
}

\details{
  Uses MQMaugment and scanMQM to scan for QTL's for each phenotype in the cross object. It is recomended to also install the snow library.
  The snow library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}

\value{
  Returns a MQMmulti object. This object is a list of scanone objects that can be plotted using plot.scanone(result[[trait]])
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu}  }

\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{ 
	\itemize{
	\item \code{\link{multitrait}} - Dataset used in the example analysis
	\item \code{\link{MQMaugment}} - Augmentationroutine to estimate missing data
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	}
}

\examples{
	#Doing a multitrait analysis
	library(MQMpackage)
	data(multitrait)
	cof <- MQMCofactorsEach(multitrait,3)
	result <- scanMQMall(multitrait,cofactors=cof,doLOG=1,b_size=3)

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
