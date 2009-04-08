\name{scanMQM}
\alias{scanMQM}

\title{ 
	scanMQM - Multiple qtl mapping for inbred crosses 
}
\description{
	Main multiQTL scanning function for the R/qtl package (for use with F2 / BC / RIL populations)
}
\usage{
	scanMQM(cross = NULL, cofactors = NULL, doLOG = 0, pheno.col = 1, REMLorML = 0, alfa = 0.02, em.iter = 1000, windowsize = 25, step.size = 5, step.min = -20, step.max = 220, n.run = 0, est.map=0,forceRIL=0,file = "MQM_output.txt",plot=TRUE)
}
\arguments{
  \item{cross}{ Cross object, an object of class \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{pheno.col}{ Phenotype under analysis is specified by the pheno.col parameter }    
  \item{doLOG}{ Indicate if the algorithm should LOG transform traitdata (DEFAULT 0 = off | 1 = on) }    
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ Windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{est.map}{ Reestimation of the marker locations on the genetic map before mapping QTLs (no map is returned the old map remains in the cross object }
  \item{dominance}{ Should dominance be taken into account ( 0= Additive  (DEFAULT)/ 1= Dominance) }  
  \item{forceRIL}{ Should we force the cross as a RIL, this option removed AB genotypes and replaces them with BB ? ( 0= NO  (DEFAULT)/ 1= YES) }  
  \item{file}{ File to which the raw results are saved (MQM_output.txt) }
  \item{plot}{ Plotting of the results ? ( T= TRUE  (DEFAULT)/ F= FALSE) }
}
\details{
Normal layout of MQM-QTL analysis
\itemize{
\item Load in experiment using formats supported by R/QTL
\item Set trait under consideration as first trait in R/QTL format
\item Fill in missing data using either the fill.geno or MQMaugmentdata routine
\item Call scanMQM to do a genomescan
\item Call scanMQM setting the n.run to get an estimate of ‘significant’ lod scores
\item Set cofactors at ‘Top scoring’ and not so ‘Top scoring’ markers
\item Rerun analysis with cofactors
\item Bootstrap analysis with cofactors to get estimates of QTL significance
}
Using ML or Reduced ML the algorithm employs a backward selection strategy to identify QTL hotspot. The algorithm passes through the following stages:
\itemize{
\item Calculating relative marker positions
\item (Re)Estimating recombinant frequencies
\item LOGlikelyhood estimations of the full model using all cofactors
\item IF(Cofactors) THEN Backward elimination of cofactors & mapQTL
\item IF(No cofactors) THEN mapQTL function
}
Afterwards the Interval mapping of the QTL model is returned to the user
}
\value{
Returns a matrix with QTLscores,Informationcontent,QTL*Informationcontent for each location on the chromosome, which can be plotted using the plot.scanone function of R/QTL or plot.MQMone to also visualize the QTL*informationcontent
The phenotype under consideration is specified by the pheno.col argument of the function.
}
\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{
\itemize{
	\item \code{\link{MQMaugment}} - Augmentationroutine to estimate missing data
	\item \code{\link{MQMCofactors}} - Learn how to set Cofactors for the MQM algorithm
	\item \code{\link{scanMQMall}} - Function called to do single trait analysis
	\item \code{\link{bootstrapMQM}} - Function called to estimate confidence interfalls for MQM
	\item \link[qtl]{scanone} - R/QTL main scanning function
	\item \link[qtl]{cim} - Composite interval mapping in R/QTL
}
}
\examples{
#Simulated F2 Population
library(MQMpackage)
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
f2result <- scanMQM(f2cross)				# Do a MQM scan of the genome
lodint(f2result,chr=3)					# LOD interval of the QTL on chr 3

#Simulated BC Population
library(MQMpackage)
bcqtl <- c(3,15,2)					# QTL at chromosome 3
data(map10)						# Mouse genome
bccross <- sim.cross(map10,bcqtl,n=100,type="bc")	# Simulate a BC Cross
bcresult <- scanMQM(bccross)				# Do a MQM scan of the genome

}
\references{
	\itemize{
		\item Jansen RC, (2007), Chapter 18 - Quantitative trait loci in inbred lines. \emph{Handbook of Stat. Genetics} 3th edition,(c) 2007 John Wiley & Sons, Ltd.
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
