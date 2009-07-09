\name{MQMfind.marker}
\alias{MQMfind.marker}

\title{ Find the significant markers after bootstrapping/permutation }

\description{
 Find the genetic markers significant after bootstrapping/permutation of the data
}

\usage{
MQMfind.marker(cross, scanMQM = NULL, perm = NULL, alpha = 0.05, verbose=FALSE)
}

\arguments{
  \item{cross}{ Cross object, An object of class \code{cross} from R/QTL with multiple phenotypes suitable for scanone and scanMQM from R/QTL. See \code{\link[qtl]{read.cross}} for details.  }
  \item{scanMQM}{ Results from either \code{\link[qtl]{scanone}} or \code{\link[qtl]{scanMQM}} }
  \item{perm}{ a \code{scanoneperm} object }
  \item{alpha}{ Threshold value, everything with significance < alpha is reported }
  \item{verbose}{If TRUE, print tracing information.}
}

%\details{
%
%}

\value{
	return a matrix with at each row a significant marker (determined from the \code{scanoneperm} object) and with columns: markername, chr and pos (Cm)
}

%\references{ 
%
%}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} }

\seealso{ 
\itemize{
	\item \code{\link[qtl]{MQM}} MQM description
	\item \code{\link[qtl]{scanMQM}} - Function called to do single trait analysis
	\item \code{\link[qtl]{bootstrap}} - Function called to estimate confidence interfalls for scanone,cim or mqm
	\item \code{\link[qtl]{MQMpermObject}} - Function called to convert MQMbootstrap into an scanoneperm object
}
}

\examples{
data(multitrait)				#Use the multitrait dataset
cof <- MQMCofactorsEach(multitrait,3)		#Set cofactors at each 3th marker
multitrait <- fill.geno(multitrait)
\dontrun{result <- bootstrap(multitrait,Funktie=mqm,cofactors=cof,pheno.col=7,n.run=50,b.size=10,doLOG=1)	# Bootstrap 50 runs in batches of 10
}\dontshow{result <- bootstrap(multitrait,Funktie=mqm,cofactors=cof,pheno.col=7,n.run=2,b.size=2,doLOG=1)}
f2perm <- MQMpermObject(result)		        #Create a permutation object
summary(f2perm)				        #What LOD score is considered significant ?
marker <- MQMfind.marker(multitrait,result[[1]],f2perm)		#Find markers with a significant QTL effect (First run is original phenotype data)
marker						#Print it to the screen
}

\keyword{ models }
