\name{mqmfind.marker}
\alias{mqmfind.marker}

\title{ Find the significant markers after bootstrapping/permutation }

\description{
 Find the genetic markers significant after bootstrapping/permutation of the data
}

\usage{
mqmfind.marker(cross, mqmscan = NULL, perm = NULL, alpha = 0.05, verbose=FALSE)
}

\arguments{
  \item{cross}{ Cross object, An object of class \code{cross} from R/QTL with multiple phenotypes suitable for scanone and mqmscan from R/QTL. See \code{\link[qtl]{read.cross}} for details.  }
  \item{mqmscan}{ Results from either \code{\link[qtl]{scanone}} or \code{\link[qtl]{mqmscan}} }
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
	\item \code{\link[qtl]{mqmscan}} - Function called to do single trait analysis
	\item \code{\link[qtl]{bootstrap}} - Function called to estimate confidence interfalls for scanone,cim or mqm
	\item \code{\link[qtl]{mqmpermObject}} - Function called to convert mqmbootstrap into an scanoneperm object
}
}

\examples{
data(multitrait)				#Use the multitrait dataset
cof <- mqmcofactorsEach(multitrait,3)		#Set cofactors at each 3th marker
multitrait <- fill.geno(multitrait)
\dontrun{result <- bootstrap(multitrait,Funktie=mqm,cofactors=cof,pheno.col=7,n.run=50,b.size=10,doLOG=1)	# Bootstrap 50 runs in batches of 10
}\dontshow{result <- bootstrap(multitrait,Funktie=mqm,cofactors=cof,pheno.col=7,n.run=2,b.size=2,doLOG=1)}
f2perm <- mqmpermObject(result)		        #Create a permutation object
summary(f2perm)				        #What LOD score is considered significant ?
marker <- mqmfind.marker(multitrait,result[[1]],f2perm)		#Find markers with a significant QTL effect (First run is original phenotype data)
marker						#Print it to the screen
}

\keyword{ models }
