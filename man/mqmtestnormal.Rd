\name{mqmtestnormal}
\alias{mqmtestnormal}

\title{ pearsons normality test used in mqm }

\description{
  Wraps a pearsons normality test from the nortest package.
}

\usage{
mqmtestnormal(cross, pheno.col = 1){
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.}
  \item{pheno.col}{ the phenotype to test for normality  }
}

%\details{
%Not available yet
%}

\value{
	Boolean indicating normality of the trait in pheno.col
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} }

\note{ 
	If you find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{
	\itemize{
	\item \code{\link[qtl]{mqmscan}} - Function called to do single trait analysis
	\item \code{\link[qtl]{mqmall}} - Function called to analyse multiple traits in the cross object
	}
}

\examples{

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

\keyword{ utilities }
