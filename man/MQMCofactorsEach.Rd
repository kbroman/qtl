\name{MQMCofactorsEach}
\alias{MQMCofactorsEach}

\title{ 
	MQMCofactorsEach - Helper function to create cofactors to be used with scanMQM 
}
\description{
  Return a cofactorslist formatted to be used with the scanMQM function
}
\usage{
MQMCofactorsEach(cross = NULL, each = 3)
}
\arguments{
  \item{cross}{ Cross object, an object of class \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
  \item{each}{ Every 'each' marker will be used as a cofactor }
}
\details{
	Not available yet
}
\value{
	An formatted list of cofactors to be used with the scanMQM algorithm
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} }

\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
}
\examples{
library(MQMpackage)
qtl <- c(3,15,3,7)						# QTL at chromosome 3
data(map10)							# Mouse genome
cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
a <- MQMCofactorsEach(cross,3)					# Set Cofactors each third marker
result <- scanMQM(cross,a)					# Do a backward model selection

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
