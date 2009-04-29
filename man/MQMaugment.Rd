\name{MQMaugment}
\alias{MQMaugment}

\title{Missing data augmentation routine for MQM}

\description{
  Estimates missing data; unlikely genotypes are kept when more likely than 1:Neglect
}

\usage{
MQMaugment(cross, pheno.col=1, maxaug=10000, maxiaug=1000, neglect=1000, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.} 
  \item{pheno.col}{Column number in the phenotype matrix which should be
   used as the phenotype.  This must be a single numeric index or a
  single character string matching one of the phenotype names.}
  \item{maxaug}{ Maximum number of (total) individuals after augmentation }
  \item{maxiaug}{ Maximum augmentations per individual }
  \item{neglect}{ Genotypes considered have a higher probability of occuring then one in neglect }
  \item{verbose}{If TRUE, tracing information is printed.}
}

\details{
Special care has 2 be taken when setting unrealistic values.
When the algorithm ends with an "DATASET TOO LARGE after augmentation",
the cross object could be corrupted.
}

\value{
  The input cross object, though with only the one phenotype indicated
  via \code{pheno.col}, and some individuals may be omitted and many
  will be repeated multiple times.  An added component \code{extra}
  contains information on which individuals are retained and repeated.
}  

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl}}

\note{
If you find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}

\seealso{
	\itemize{
	\item \code{\link[qtl]{scanMQM}} - Function called to do single trait analysis
	\item \code{\link[qtl]{scanall}} - Function called to analytse multiple traits in the cross object
	}
}

\examples{
qtl <- c(3,15,3,7)							# QTL at chromosome 3
data(map10)								# Mouse genome
cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)			# Simulate a Cross
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
