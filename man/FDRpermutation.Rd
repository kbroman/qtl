\name{FDRpermutation}
\alias{FDRpermutation}

\title{ Estimates of FDR when doing multiple trait QTL analysis }

\description{
   Using a randomization approaches to obtain estimates of FDR when doing multiple trait QTL analysis
}

\usage{
	FDRpermutation(cross, Funktie=scanall, thresholds=c(1,2,3,4,5,7,10,15,20), n.perm=10, verbose=TRUE, \dots)
}

\arguments{
	\item{cross}{ Cross object, An object of class \code{cross} from R/QTL with multiple phenotypes suitable for scanone and mqmscan from R/QTL. See \code{\link[qtl]{read.cross}} for details. }
	\item{Funktie}{ Function to use when mappingQTL's (either scanall,cimall or mqmall)}
	\item{thresholds}{ False discovery rate is calculated for peaks above these LOD thresholds (DEFAULT=Range from 1 to 20, using 10 thresholds) }
	\item{n.perm}{ Number of permutation runs (DEFAULT=10, should be 1000 to 5000) }
	\item{verbose}{If TRUE, print tracing information}
	\item{\dots}{Parameters passed through to the
	  \code{\link[qtl]{scanall}}, \code{\link[qtl]{cimall}} or
	    \code{\link[qtl]{mqmall}} functions} 
}

\details{
  Analysis of \code{\link[qtl]{scanall}}, \code{\link[qtl]{cimall}} or
	    \code{\link[qtl]{mqmall}} to scan for QTL in shuffled/randomized data. It is recommended to also install the \code{snow} library.
  The \code{snow} library allows calculations to run on multiple cores or even scale it up to an entire cluster, thus speeding up calculation.
}

\value{
	Returns a matrix with FP FN and FDR for the selected thresholds
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} ; Karl W. Broman \email{kbroman@biostat.wisc.edu} }

\note{ 
	If you find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{
	\itemize{
	\item \code{\link[qtl]{MQMCofactors}} - Learn how to set Cofactors for the MQM algorithm
	\item \code{\link[qtl]{scanall}} - Function called to analyse multiple traits in the cross object (
	\item \code{\link[qtl]{cimall}} - Function called to analyse multiple traits in the cross object (Using CIM)
	\item \code{\link[qtl]{mqmall}} - Function called to analyse multiple traits in the cross object (Using MQM)
	}
}

\examples{
data(multitrait)				#Use the multitrait dataset
\dontshow{multitrait <- subset(multitrait, chr=1, ind=!apply(multitrait$pheno, 1, function(a) any(is.na(a))))}
multitrait <- calc.genoprob(multitrait)
result <- FDRpermutation(multitrait, threshold=10.0, n.perm=5, method="hk")	#Calculate the thresholds
}

\references{
	\itemize{
		\item Jansen RC, (2007) Quantitative trait loci in inbred lines. Chapter 18 of \emph{Handbook of Stat. Genetics} 3rd edition. John Wiley & Sons, Ltd.
		\item Tierney, L., Rossini, A., Li, N., and Sevcikova, H. (2004), The snow Package: Simple Network of Workstations. Version 0.2-1. 	
		\item Rossini, A., Tierney, L., and Li, N. (2003), Simple parallel statistical computing.  \emph{R. UW Biostatistics working paper series} University of Washington. \bold{193}
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388--391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 1447--1455.
		\item Churchill, G. A. and Doerge, R. W. (1994) Empirical threshold values for quantitative trait mapping. \emph{Genetics} \bold{138}, 963--971.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205--211.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 1--38.
		\item Bruno M. Tesson1, Ritsert C. Jansen (2009) Chapter 3.7. Determining the significance threshold \emph{eQTL Analysis in Mice and Rats} \bold{1}, 20--25
 	}
}

\keyword{ models }
