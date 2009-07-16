\name{MQMCofactors}
\alias{MQMCofactors}

\title{Create cofactors to be used with mqmscan }

\description{
  Return a cofactorslist formatted to be used with the mqmscan function
}

\usage{
MQMCofactors(cross, cofactors, sexfactors, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
  \item{cofactors}{ List of integer describing which marker should be treated as a cofactor }
  \item{sexfactors}{ List of integer describing which marker should be treated as a sexfactor }
  \item{verbose}{If TRUE, print tracing information.}
}

%\details{
%Not available yet
%}

\value{
An formatted list of cofactors to be used with the \code{\link[qtl]{mqmscan}} algorithm
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} }

\note{ 
	If you find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{
	\itemize{
	\item \code{\link[qtl]{mqmscan}} - Function called to do single trait analysis
	\item \code{\link[qtl]{scanall}} - Function called to analytse multiple traits in the cross object
	}
}

\examples{
qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
a <- MQMCofactors(cross,c(20,30,80,90),c(186,187))	# Set Cofactors on marker 20,30,80 & 90 and sexfactors on 186,187
result <- mqmscan(cross,a)				# Do a backward model selection
}

\references{
	\itemize{
		\item Jansen RC, (2007) Quantitative trait loci in inbred lines. Chapter 18 of \emph{Handbook of Stat. Genetics} 3rd edition. John Wiley & Sons, Ltd.
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388--391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 1447--1455.
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963--971.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205--211.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 1--38.
	}
}

\keyword{ models }
