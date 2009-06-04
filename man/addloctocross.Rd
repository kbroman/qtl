\name{addloctocross}
\alias{addloctocross}

\title{Adds phenotypes location into a cross object}

\description{
	Adds genomic locations of phenotype (phenotypes that have genomic locations usually come from eQTL or pQTL studies) into a cross object
}

\usage{
	addloctocross(cross,locations=NULL, locfile="locations.txt")
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.} 
  \item{locations}{Rvariable holding location information about the traits.}
  \item{locfile}{ if there a no \code{locations} specified data is loaded from a file specified by this parameter, see the details section for the layout of the file.}
}

\details{
inputfile layout:
Num	Name	Chr	cM
1	X3.Hydroxypropyl	4	50.0
Num is the number of the phenotype in the cross object
Name is the name of the phenotype (will be checked against the name already in the cross object at position num
Chr Chromosome
cM	postion from start of chromosome in cM
}

\value{
  The input cross object, with the locations added as an aditional component \code{locations}
}  

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl}}

\note{
If you find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}

\seealso{
	\itemize{
		\item \code{\link[qtl]{CisTransPlot}} - Function called to analyze multiple traits in the cross object
	}
}

\examples{
% See the help file for CisTransPlot
\dontrun{
	data(multitrait)
	data(locations)
	multiloc <- addloctocross(multitrait,locations)
	results <- scanall(multiloc)
	CisTransPlot(results, multiloc, 5, FALSE, TRUE)
}
}

%\references{
%
%}

\keyword{ utilities }
