\name{bootstrapMQM}
\alias{bootstrapMQM}
\alias{permutateMQM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ bootstrapMQM - Bootstrapping approach to obtain estimates of QTL LOD score significance }
\description{
   Bootstrapping approach to obtain estimates of QTL significance
}
\usage{
bootstrapMQM(cross = NULL, cofactors = NULL, Phenot = 1, REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.run = 10, file = "MQM_output.txt", doLOG = 0, est.map = 0, dominance = 0, n.clusters = 2, parametric = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object, An object of class \code{cross} from R/QTL with multiple phenotypes suitable for scanone and scanMQM from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{n.run}{ Number of bootstrapping runs (DEFAULT=10) }
  \item{parametric}{ What kind of bootstrapping should be done ( 0= Random Shuffle (DEFAULT)/ 1= Parametric bootstrapping ) }  
  \item{Phenot}{ Phenotype to do thepull from the cross object }
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or Maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ Windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{doLOG}{ Indicate if the algorithm can LOG transform traits with variation > 10.000 (DEFAULT 0 = off | 1 = on) }
  \item{est.map}{ Reestimation of the genetic map (est.rf in R/QTL) }
  \item{dominance}{ Should dominance be taken into account ( 0= Additive  (DEFAULT)/ 1= Dominance) }
  \item{n.clusters}{ Number of child processes to split the job into }
  \item{b_size}{ Batch size. The entire set is split in jobs. Each job contains b_size # of traits per job }
  
}
\details{
  scanMQM to scan for QTL's in bootstrapped data. It is recomended to also install the snow library.
  The snow library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}
\value{
	Returns a MQMmulti object. this object is a list of scanone objects that can be plotted using plot.scanone(result[[trait]])
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu} }

\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{
	\itemize{
	\item \code{\link{MQMCofactors}} - Learn how to set Cofactors for the MQM algorithm
	\item \code{\link{scanMQMall}} - Function called to do single trait analysis
	}
}
\examples{
library(MQMpackage)
data(multitrait)				#Use the multitrait dataset
cof <- MQMCofactorsEach(multitrait,3)		#Set cofactors at each 3th marker
f2result <- bootstrapMQM(multitrait,cofactors=cof,pheno.col=7,n.run=50,b_size=20,doLOG=1)	# Bootstrap 200 runs in batches of 20

#EXTRA
f2perm <- MQMpermObject(f2result)			#Create a permutation object
summary(f2perm)						#Get Significant LOD thresholds

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
