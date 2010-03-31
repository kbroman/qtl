\name{mqmfind.marker}
\alias{mqmfind.marker}

\title{ Fetch significant markers after permutation analysis }

\description{
 Fetch significant makers after permutation analysis.  These markers can be
 used as cofactors for model selection in a forward stepwise approach.
}

\usage{
mqmfind.marker(cross, mqmscan = NULL, perm = NULL, alpha = 0.05, verbose=FALSE)
}

\arguments{
  \item{cross}{ 
An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
  \item{mqmscan}{ Results from either \code{\link{scanone}} or \code{\link{mqmscan}} }
  \item{perm}{ a \code{scanoneperm} object }
  \item{alpha}{ Threshold value, everything with significance < alpha is reported }
  \item{verbose}{
Display more output on verbose=TRUE % \verbose
  }
}

\value{
  returns a matrix with at each row a significant marker (determined from the
  \code{scanoneperm} object) and with columns: markername, chr and pos (Cm) 
}

\author{ 
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{ 
\itemize{
	\item \code{\link{mqmprocesspermutation}} - Function called to convert results from an mqmpermutation into an scanoneperm object
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
cof <- mqmsetcofactors(multitrait,3)		#Set cofactors at each 3th marker
multitrait <- fill.geno(multitrait)
\dontrun{result <- bootstrap(multitrait,scanfunction=mqmscan,cofactors=cof,pheno.col=7,n.run=50,b.size=10,logtransform=TRUE)	# Bootstrap 50 runs in batches of 10
}\dontshow{result <- bootstrap(multitrait,scanfunction=mqmscan,cofactors=cof,pheno.col=7,n.run=2,b.size=2,logtransform=TRUE)}
f2perm <- mqmprocesspermutation(result)		        #Create a permutation object
summary(f2perm)				        #What LOD score is considered significant ?
marker <- mqmfind.marker(multitrait,result[[1]],f2perm)		#Find markers with a significant QTL effect (First run is original phenotype data)
marker						#Print it to the screen
}

\keyword{ models }
