\name{mqmall}
\alias{mqmall}

\title{QTL analysis using mqm on multiple phenotypes in a cross object}

\description{
	QTL analysis using mqm on multiple phenotypes in a cross object uses the SNOW library for faster computations
}

\usage{
	mqmall(cross, multiC=TRUE, n.clusters = 1,b.size=10, FF=0,\dots, plot=FALSE, verbose=FALSE)
}

\arguments{
  \item{cross}{ Cross object, an object of class 
    \code{cross} from R/QTL. See \code{\link[qtl]{read.cross}} for details. 
  }
  \item{multiC}{ Use multiple cores (only if the package SNOW is available)
  }
  \item{n.clusters}{ Number of child processes to split the job into. with multiC 
  is TRUE this determins how many cores are used at a given time. For optimal performance use
  2 on a duocore, 3 on a quadcore and 5 on a hexacore.
  }
  \item{b.size}{ Batch size. The entire set is split in 
    jobs. Each job contains b.size number of traits per job 
  }
  \item{\dots}{ Parameters passed through to the \code{\link[qtl]{mqm}} function used in 
    scanning for QTLs
  }
}

\details{
  Uses the mqm function to scan for QTL's for each 
  phenotype in the cross object. It is recomended 
  to also install the snow library. The snow library 
  allows calculations to run on multiple cores or 
  even scale it up to an entire cluster, thus speeding up calculation.
}

\value{
  Returns a MQMmulti object. This object is a list of \code{\link[qtl]{scanone}}
  objects that can be plotted using \code{plot.scanone(result[[trait]])} or using \code{mqmplot_multitrait(result)}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} ; Karl W. Broman \email{kbroman@biostat.wisc.edu}  }

\note{ 
	If you find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{ 
	\itemize{
	\item \code{\link[qtl]{multitrait}}  Dataset used in the example analysis
	\item \code{\link[qtl]{fill.geno}}  Augmentationroutine to estimate missing data	
	\item \code{\link[qtl]{mqmaugment}}  Augmentationroutine to estimate missing data
	\item \code{\link[qtl]{mqm}}  Function called to do single trait QTL analysis
	}
}

\examples{
	#Doing a multitrait analysis
	data(multitrait)
  \dontshow{multitrait <- subset(multitrait, chr=1:2, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))
  multitrait$pheno <- multitrait$pheno[,1:3]}
  multitrait <- calc.genoprob(multitrait)
	cof <- mqmsetcofactors(multitrait,3)
  multitrait <- fill.geno(multitrait)
	result <- mqmall(multitrait,cofactors=cof,b.size=5)
}

\references{
	\itemize{
		\item Jansen RC, (2007) Quantitative trait loci in inbred lines. Chapter 18 of \emph{Handbook of Stat. Genetics} 3rd edition. John Wiley & Sons, Ltd.
		\item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.21. 	
		\item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 14471455.
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963971.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205211.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 138.
	}
}

\keyword{ models }