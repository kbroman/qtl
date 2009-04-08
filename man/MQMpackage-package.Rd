\name{MQMpackage-package}
\alias{MQMpackage}
\alias{MQM}
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
	Package: \tab MQMpackage\cr
	Version: \tab 1.0\cr
	Date: \tab 30 January 2009\cr
	Depends: \tab R (>= 1.8.0), qtl\cr
	License: \tab GPL version 2 or later\cr
	LazyLoad: \tab yes\cr
	LazyData: \tab yes\cr
	Packaged: \tab Mon Feb  9 13:36:40 2009; F101715\cr
	Built: \tab R 2.8.1; i386-pc-mingw32; 2009-02-09 13:36:42; windows\cr
	}

	Index:
	\preformatted{
	scanMQM                 scanMQM - Multiple QTL mapping for F2 crosses
	scanMQMall              scanMQMall - Multitrait analysis of a cross object
	bootstrapMQM		bootstrapMQM - Bootstrapping approach to obtain estimates of QTL significance
	MQMaugment              MQMaugment - Missing data augmentation routine for MQM
	MQMCofactors            MQMCofactors - Helper function to create cofactors to be used with scanMQM
	MQMCofactorsEach        MQMCofactorsEach - Helper function to create cofactors to be used with scanMQM
	CrossFromMolgenis       CrossFromMolgenis - Helperfunction to create a cross-object from a molgenis database
	ResultsFromMolgenis     ResultsFromMolgenis - Helperfunction to retrieve results from a molgenis database
	ResultsToMolgenis       ResultsToMolgenis - Helperfunction to store results into a molgenis database
	multitrait              Crossobject from R/QTL with multiple traits
	PipelineMolgenis		Automatic MQM-QTL analysis for Molgenis database systems	
	}
	Normal layout of MQM-QTL analysis
	\itemize{
		\item Load in experiment using formats supported by R/QTL. See \code{\link[qtl]{read.cross}} for details.
		\item Set trait under consideration as first trait in R/QTL format
		\item Fill in missing data using either the fill.geno or MQMaugmentdata routine
		\item Call scanMQM to do a genomescan
		\item Call scanMQM setting the n.run to get an estimate of ‘significant’ lod scores
		\item Set cofactors at ‘Top scoring’ and not so ‘Top scoring’ markers
		\item Rerun analysis with cofactors
		\item Bootstrap analysis with cofactors to get estimates of QTL significance
	}
	The Molgenispipeline is comprised of two main parts:
	\itemize{
		\item 1A: Data retrieval from a molgenis database 
		\item 1B: Formatting into a \code{cross} object. See \code{\link[qtl]{read.cross}} for details. 
		\item 2A.1: MQM QTL Model selection, with specified cofactors
		\item 2A.2: MQM QTL mapping using the selected model
		\item 2B: Storage of QTLmap into the database
	}
}
\author{
Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} ; Karl W. Brohamn \email{kbroman@biostat.wisc.edu} ; Rutger Brouwer \email{R.W.W.Brouwer@rug.nl} ; M.A. Swertz \email{m.a.swertz@rug.nl}

Maintainer: Danny Arends \email{Danny.Arends@gmail.com}
}
\references{
	\itemize{
		\item Jansen RC, (2007), Chapter 18 - Quantitative trait loci in inbred lines. \emph{Handbook of Stat. Genetics} 3th edition,(c) 2007 John Wiley & Sons, Ltd.
		\item Swertz MA, Jansen RC. (2007), Beyond standardization: dynamic software infrastructures for systems biology. \emph{Nat Rev Genet.} \bold{3}, 235--243.
		\item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.2-1. 	
		\item Swertz MA, De Brock EO, Van Hijum SA, De Jong A, Buist G, Baerends RJ, Kok J, Kuipers OP, Jansen RC. (2004), Molecular Genetics Information System (MOLGENIS): alternatives in developing local experimental genomics databases. \emph{Bioinformatics},\bold{13}, 2075--2083.
		\item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388-391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 1447-1455.
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963--971.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205--211.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 1--38.
	}
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }

\seealso{
	\itemize{
	\item	\code{\link{MQMaugment}} - Dataaugmentation routine to estimate missing data
	\item	\code{\link{scanMQM}} - Function called to do single trait analysis
	\item	\code{\link{scanMQMall}} - Function called to do multiple trait analysis
	\item \code{\link{PipelineMolgenis}} - Automatic MQM-QTL analysis for Molgenis database systems		
	}
}
\examples{
	library(MQMpackage)
	data(multitrait)
	mt_est <- MQMaugment(multitrait)
	result <- scanMQM(mt_est)
	plot(result)
}
