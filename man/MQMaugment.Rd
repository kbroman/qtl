\name{MQMaugment}
\alias{MQMaugment}

\title{MQM data augmentation}

\description{
  Estimates missing data. Unlikely genotypes are kept when more likely than the range 1:Neglect
}

\usage{
MQMaugment(cross, pheno.col=1, maxaug=1000, maxiaug=10, neglect=10,
           verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}}.} 
  \item{pheno.col}{Column number in the phenotype matrix. This can be an index
	number or a character string matching one of the phenotype names. If pheno.col is
	an array the first item is used}
  \item{maxaug}{ Maximum number of (expanded) individuals after augmentation }
  \item{maxiaug}{ Maximum augmentations per individual }
  \item{neglect}{ Ascertain genotypes below this probability are ignored }
  \item{verbose}{If TRUE, tracing information is printed.}
}

\details{
}

\value{
	The cross object the single phenotype \code{pheno.col}, without dropped
	individuals and with augmented individuals (many will be repeated multiple
	times). An added component \code{extra} contains information on which
	individuals are retained and repeated.  
}  

\author{R.C. Jansen ; Danny Arends}

\note{
The sex chromosome 'X' is dropped before augmentation.

With unrealistic values for augmentation a warning is shown and the cross
object may have been corrupted (FIXME).  }

\seealso{
	\itemize{
	\item \code{\link[qtl]{scanMQM}} - Function called to do single trait analysis
	\item \code{\link[qtl]{scanall}} - Function called to analytse multiple traits in the cross object
	}
}

\examples{
qtl <- c(3,15,3,7)							# QTL at chromosome 3
data(map10)								# Mouse genome
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01)             # Simulate a Cross (just autosomes)
crossaug <- MQMaugment(cross)						# Augmentation of the data
result <- scanMQM(crossaug)						# Do a MQM scan of the genome
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
