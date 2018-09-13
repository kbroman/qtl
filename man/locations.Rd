\name{locations}
\alias{locations}

\docType{data}

\title{
    Genetic locations of traits for the multitrait dataset
}

\description{
    A table with genetic locations of the traits in the \code{\link{multitrait}} dataset
}
\usage{data(locations)}
\format{
    Each row is a trait with the following information:
    \code{Name},	Name of the trait (will be checked against the name in the cross object
    \code{Chr},		Chromosome of the trait
    \code{cM},		Location in cM from the start of the chromosome
}
%\details{
%
%}
\source{
  Additional information from the Arabidopsis RIL selfing experiment with Landsberg erecta (Ler) and Cape Verde Islands (Cvi) with 162 individuals scored (with errors at) 117 markers.
  Dataset obtained from GBIC - Groningen BioInformatics Centre
}
\references{
\itemize{

 \item Keurentijes JJB, Fu J, de Vos CHR,Lommen A, Jansen RC et al
 (2006), The genetics of plant metabolism. \emph{Nature Genetics}
 \bold{38}, 842--849.

 \item Alonso-Blanco C., Peeters, A. J. and Koornneef, M. (2006)
 Development of an AFLP based linkage map of Ler, Col and Cvi
 Arabidopsis thaliana ecotypes and construction of a Ler/Cvi recombinant
 inbred line population. \emph{Plant J.} \bold{14}(2), 259--271.
}

}

\examples{
% See the help file for mqmplot.cistrans
\dontrun{
    data(multitrait)
    data(locations)
    multiloc <- addloctocross(multitrait,locations)
    results <- scanall(multiloc)
    mqmplot.cistrans(results,multiloc, 5, FALSE, TRUE)
}
}

\seealso{\code{\link{multitrait}} }

\keyword{datasets}
