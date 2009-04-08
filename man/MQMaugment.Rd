\name{MQMaugment}
\alias{MQMaugment}

\title{ 
	MQMaugment - Missing data augmentation routine for MQM 
}
\description{
  Estimates missing data, unlikely genotypes are kept when more likely than 1:Neglect
}
\usage{
MQMaugment(cross = NULL, maxaug = 10000, maxiaug = 1000, neglect = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{Phenot}{ Phenotype to pull from the cross object }  
  \item{maxaug}{ Maximum number of (total) individuals after augmentation }
  \item{maxiaug}{ Maximum augmentations per individual }
  \item{neglect}{ Genotypes considered have a higher probability of occuring then one in neglect }
}
\details{
Special care has 2 be taken when setting unrealistic values.
When the algorithm ends with an "DATASET TOO LARGE after augmentation", the crossfile could be corrupted.
}
\value{
Cross object, an object of class \code{cross} from R/QTL with 1 phenotype (indicated by setting the Phenot parameter) suitable for scanMQM from R/QTL. See \code{\link[qtl]{read.cross}} for details.
}
\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl}}
\note{
If u find anything wrong with the algorithm please notify me at: Danny.Arends@gmail.com
}
\seealso{
	\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do multiple trait analysis
	}
}
\examples{
	library(MQMpackage)
	qtl <- c(3,15,3,7)							# QTL at chromosome 3
	data(map10)								# Mouse genome
	cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)			# Simulate a Cross
	cross_good <- MQMaugment(cross)						# Augmentation of the data
	result <- scanMQM(cross_good)						# Do a MQM scan of the genome
	
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
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
