\name{MQMpackage-package}
\alias{MQMpackage}
\docType{package}
\title{
	Tool for analyzing QTLs in experimental crosses
}
\description{
	Multi-trait MarkerQTLMarker(MQM) Analysis of experimental crosses to identify QTLs contributing to variation in quantitative traits.
	Currently supported crosses:
	\itemize{
		\item F2 Intercross
		\item Backcross
		\item Recombinant inbred line (Sibling mating) 
		\item Recombinant inbred line (Self-fertilization)
	}
}
\details{
	\tabular{ll}{
	Package: \tab MQMpackage, subpackage of the R/QTL package\cr
	Date: \tab 2008-2009\cr
	Depends: \tab R (>= 1.8.0), qtl\cr
	License: \tab GPL version 2 or later\cr
	LazyLoad: \tab yes\cr
	LazyData: \tab yes\cr
	Packaged: \tab Mon Feb  9 13:36:40 2009; F101715\cr
	Built: \tab R 2.8.1; i386-pc-mingw32; 2009-02-09 13:36:42; windows\cr
	}

	Index:
	\preformatted{
	mqmscan                 mqmscan - Multiple QTL mapping for F2 crosses
	scanall					scanall - Function called to analytse multiple traits in the cross object
	bootstrap			bootstrap - Bootstrapping approach to obtain estimates of QTL significance
	MQMaugment              MQMaugment - Missing data augmentation routine for MQM
	MQMCofactors            MQMCofactors - Helper function to create cofactors to be used with mqmscan
	MQMCofactorsEach        MQMCofactorsEach - Helper function to create cofactors to be used with mqmscan
	multitrait              Crossobject from R/QTL with multiple traits
	}
	Normal layout of MQM-QTL analysis
	\itemize{
		\item Load in experiment using formats supported by R/QTL. See \code{\link[qtl]{read.cross}} for details.
		\item Set trait under consideration as first trait in R/QTL format
		\item Fill in missing data using either the fill.geno or MQMaugmentdata routine
		\item Call mqmscan to do a genomescan
		\item Call mqmscan setting the n.run to get an estimate of "significant" lod scores
		\item Set cofactors at "Top scoring" and not so "Top scoring" markers
		\item Rerun analysis with cofactors
		\item Bootstrap analysis with cofactors to get estimates of QTL significance
	}
}
\author{
Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} ; Karl W. Broman \email{kbroman@biostat.wisc.edu} ; Rutger Brouwer \email{R.W.W.Brouwer@rug.nl} ; M.A. Swertz \email{m.a.swertz@rug.nl}

Maintainer: Danny Arends \email{Danny.Arends@gmail.com}
}
\references{
	\itemize{
		\item Jansen RC, (2007) Quantitative trait loci in inbred lines. Chapter 18 of \emph{Handbook of Stat. Genetics} 3rd edition. John Wiley & Sons, Ltd.
		\item Swertz MA, Jansen RC. (2007), Beyond standardization: dynamic software infrastructures for systems biology. \emph{Nat Rev Genet.} \bold{3}, 235--243.
		\item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.2-1. 	
		\item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388--391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 1447--1455.
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963--971.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205--211.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 1--38.
	}
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }

\seealso{
	\itemize{
	\item	\code{\link[qtl]{MQMaugment}} - Dataaugmentation routine to estimate missing data
	\item	\code{\link[qtl]{mqmscan}} - Function called to do single trait analysis
	\item \code{\link[qtl]{scanall}} - Function called to analytse multiple traits in the cross object
	}
}
\examples{
	data(multitrait)
\dontshow{multitrait <- subset(multitrait, chr=1:2)}        
	mt.est <- MQMaugment(multitrait)
	result <- mqmscan(mt.est, doLOG=1)
	plot(result)
}
