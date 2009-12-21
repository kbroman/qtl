\name{multitrait}
\alias{multitrait}
\alias{ArabidopsisLERCVI}

\docType{data}

\title{ 
	Example Cross object from R/QTL with multiple traits.
}

\description{
	Cross object from R/QTL, an object of class \code{cross} from R/QTL. See \code{\link{read.cross}} for details.
}
\usage{data(multitrait)}
\format{
	Cross object from R/QTL
}
\details{
	Arabidopsis Recombinant Inbred Line by selfing
}
\source{
  Part of the Arabidopsis RIL selfing experiment with Landsberg erecta (Ler) and Cape Verde Islands (Cvi) with 162 individuals scored (with errors at) 117 markers.
  Dataset obtained from GBIC - Groningen BioInformatics Centre
}
\references{
\itemize{
 
 \item Keurentijes JJB, Fu J, de Vos CHR,Lommen A, Jansen RC et al
 (2006), The genetics of plant metabolism. \emph{Nature Genetics.}
 \bold{38}, 842--849. 
 
 \item Alonso-Blanco C., Peeters, A. J. and Koornneef, M. (2006)
 Development of an AFLP based linkage map of Ler, Col and Cvi
 Arabidopsis thaliana ecotypes and construction of a Ler/Cvi recombinant
 inbred line population. \emph{Plant J.} \bold{14}(2), 259--271.
}

}

\examples{
data(multitrait)						#Load dataset
multitrait <- fill.geno(multitrait)                     # impute missing genotype data
result <- scanall(multitrait,Funktie=mqm,doLOG=1)	#Analyse all 24 traits
}

\keyword{datasets}
