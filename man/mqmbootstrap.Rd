\name{bootstrap}
\alias{bootstrap}
\alias{bootstrapcim}
\alias{mqmbootstrap}
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
	bootstrap(cross, Funktie=scanone, pheno.col=1,
	          multiC=TRUE, n.run=10, b.size=10, file="MQM_output.txt",
   	          n.clusters=1, bootmethod=0, plot=FALSE, verbose=FALSE, \dots)
	bootstrapcim(\dots)
	mqmbootstrap(\dots)
}

\arguments{
	\item{cross}{ Cross object, An object of class \code{cross} from R/QTL with multiple phenotypes suitable for scanone and mqmscan from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
	\item{Funktie}{ Function to use when mappingQTL's (either scanone,cim or mqm)}
	\item{pheno.col}{ Phenotype under investigation (from the cross object) }
	\item{multiC}{ Use multicore (if available)}  
	\item{n.run}{ Number of bootstrapping runs (DEFAULT=10, should be 1000 to 5000) }
	\item{b.size}{ Batch size. The entire set is split in jobs. Each job contains b.size number of traits per job }
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
	Returns a mqmmulti object. this object is a list of scanone objects that can be plotted using \code{plot.scanone(result[[trait]])}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} ; Karl W. Broman \email{kbroman@biostat.wisc.edu} }

\seealso{
	\itemize{
	\item \code{\link[qtl]{mqmcofactors}} - Learn how to set Cofactors for the MQM algorithm
	\item \code{\link[qtl]{scanall}} - Function called to analyse multiple traits in the cross object
	\item \code{\link[qtl]{scanone}} - Function called to do single trait analysis
	\item \code{\link[qtl]{cim}} - Function called to do single trait analysis
	\item \code{\link[qtl]{mqm}} - Function called to do single trait analysis
	}
}

\examples{
data(multitrait)				#Use the multitrait dataset
\dontshow{multitrait <- subset(multitrait, chr=1:2, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
multitrait$pheno <- log(multitrait$pheno)}
multitrait <- calc.genoprob(multitrait)
\dontrun{result <- bootstrap(multitrait,pheno.col=7, b.size=5)
result <- bootstrapcim(multitrait,pheno.col=7, b.size=5)
}\dontshow{result <- bootstrap(multitrait,pheno.col=7, n.run=2, b.size=2)
result <- bootstrapcim(multitrait,pheno.col=7, n.run=2, b.size=2)}
multitrait <- fill.geno(multitrait)
cof <- mqmcofactorsEach(multitrait,3)		#Set cofactors at each 3th marker
\dontrun{f2result <- mqmbootstrap(multitrait,cofactors=cof,pheno.col=7, b.size=5,verbose=FALSE)}
\dontshow{f2result <- mqmbootstrap(multitrait,cofactors=cof,pheno.col=7, n.run=2,b.size=2,verbose=FALSE)}
#EXTRA
f2perm <- mqmpermObject(f2result)			#Create a permutation object
summary(f2perm)						#Get Significant LOD thresholds
}

\keyword{ models }
