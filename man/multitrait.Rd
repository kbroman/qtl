\name{multitrait}
\alias{multitrait}
\docType{data}

\title{ 
	Example Cross object from R/QTL with multiple traits
}

\description{
	Cross object from R/QTL, an object of class \code{cross} from R/QTL. See \code{\link{read.cross}} for details.
}
\usage{data(multitrait)}
\format{
	Cross object from R/QTL
}
\details{
	Arabidopsis recombinant inbred lines by selfing.  There are 162
	lines, 24 phenotypes, and 117 markers on 5 chromosomes.
}
\source{
  Part of the Arabidopsis RIL selfing experiment with Landsberg erecta (Ler) and Cape Verde Islands (Cvi) with 162 individuals scored (with errors at) 117 markers.
  Dataset obtained from GBIC - Groningen BioInformatics Centre
}
\references{
\itemize{
 
 \item Keurentjes, J. J. and Fu, J. and de Vos, C. H. and Lommen, A. and Hall, R. D. and Bino, R. J. and van der Plas, L. H. and Jansen, R. C. and Vreugdenhil, D. and Koornneef, M.
 (2006), The genetics of plant metabolism. \emph{Nature Genetics.}
 \bold{38}-7, 842--849. 

 \item Alonso-Blanco, C. and Peeters, A. J. and Koornneef, M. and Lister, C. and Dean, C. and van den Bosch, N. and Pot, J. and Kuiper, M. T. (1998),
 Development of an AFLP based linkage map of Ler, Col and Cvi
 Arabidopsis thaliana ecotypes and construction of a Ler/Cvi recombinant
 inbred line population}. \emph{Plant J.} \bold{14}(2), 259--271.

}

\examples{
data(multitrait)					# Load dataset
multitrait <- fill.geno(multitrait)                     # impute missing genotype data
\dontshow{multitrait <- subset(multitrait, chr=1)
multitrait$pheno <- multitrait$pheno[,1:3]}
result <- mqmscanall(multitrait, logtransform=TRUE)	# Analyse all 24 traits
}

\keyword{datasets}
