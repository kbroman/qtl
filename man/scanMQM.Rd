\name{scanMQM}
\alias{scanMQM}
\alias{mqm}

\title{Multiple qtl mapping for inbred crosses}


\description{
	Main multiQTL scanning function (for use with F2 / BC / RIL populations)
}

\usage{
	scanMQM(cross, cofactors, pheno.col = 1, REMLorML = 0, 
                alfa = 0.02, em.iter = 1000, windowsize = 25.0, step.size = 5.0, 
                step.min = -20.0, step.max = 220.0, file = "MQM_output.txt",
                doLOG=0, est.map=0, dominance=0, plot=FALSE, forceRIL=0, verbose=FALSE)
}

\arguments{
  \item{cross}{ Cross object, an object of class \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{pheno.col}{ Phenotype under analysis is specified by the pheno.col parameter }    
  \item{REMLorML}{ Lets the algorithm use reduced maximum likelyhood (0) or maximum likelyhood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ Windowsize used when mappingQTL locations }
  \item{step.size}{ Step.size used for the generation of output QTL }
  \item{step.min}{ Begin this many Cm before the chromosome }
  \item{step.max}{ Maximum length in Cm }
  \item{file}{ File to which the raw results are saved.}
  \item{doLOG}{ Indicate if the algorithm should LOG transform traitdata (DEFAULT 0 = off | 1 = on) }    
  \item{est.map}{ Reestimation of the marker locations on the genetic map before mapping QTLs (no map is returned the old map remains in the cross object }
  \item{dominance}{ Should dominance be taken into account ( 0= Additive  (DEFAULT)/ 1= Dominance) }  
  \item{forceRIL}{ Should we force the cross as a RIL, this option removed AB genotypes and replaces them with BB ? ( 0= NO  (DEFAULT)/ 1= YES) }  
  \item{plot}{ If TRUE, plot the results.}
  \item{verbose}{If TRUE, print tracing information.}
}

\details{
Normal layout of MQM-QTL analysis
\itemize{
\item Load in experiment using formats supported by R/QTL
\item Set trait under consideration as first trait in R/QTL format
\item Fill in missing data using either the fill.geno or MQMaugmentdata routine
\item Call scanMQM to do a genomescan
\item Call scanMQM setting the n.run to get an estimate of "significant" lod scores
\item Set cofactors at "Top scoring" and not so "Top scoring" markers
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

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} ; Karl W. Broman \email{kbroman@biostat.wisc.edu} }

\note{ 
	If you find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{
\itemize{
	\item \code{\link[qtl]{MQMaugment}} - Augmentationroutine to estimate missing data
	\item \code{\link[qtl]{MQMCofactors}} - Learn how to set Cofactors for the MQM algorithm
	\item \code{\link[qtl]{scanall}} - Function called to analytse multiple traits in the cross object
	\item \code{\link[qtl]{bootstrap}} - Function called to estimate confidence interfalls for MQM
	\item \code{\link[qtl]{scanone}} - R/QTL main scanning function
	\item \code{\link[qtl]{cim}} - Composite interval mapping in R/QTL
}
}

\examples{
#Simulated F2 Population
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
f2result <- scanMQM(f2cross)				# Do a MQM scan of the genome
lodint(f2result,chr=3)					# LOD interval of the QTL on chr 3

#Simulated BC Population
bcqtl <- c(3,15,2)					# QTL at chromosome 3
data(map10)						# Mouse genome
bccross <- sim.cross(map10,bcqtl,n=100,type="bc")	# Simulate a BC Cross
bcresult <- scanMQM(bccross)				# Do a MQM scan of the genome
}

\keyword{ models }
