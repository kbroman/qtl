\name{MQM}
\alias{MQM}

\title{Multiple QTL Model (MQM) mapping}

\description{
Overview of the MQM mapping functions
}

\usage{
}

\section{Introduction}{
MQM is an implementation of multiple QTL model mapping (or Marker-QTL-Marker)
as described by Ritsert Jansen (see references below). A genome scan is
performed by moving a QTL along the chromosomes while using a preselected set
of markers as cofactors. If the QTL gets close to a cofactor the cofactor gets
dropped. The backward elimination procedure starts with a multiple-regression
model on an augmented dataset, using a full set of cofactors evenly spread over
the genome. A maximum likelihood (ML) approach is used for cofactor selection.
}

\author{R.C. Jansen ; Danny Arends ; Karl Broman \email{kbroman@biostat.wisc.edu} ; Pjotr Prins}

\note{
Currently MQM only supports F2 / BC / RIL crosses. Also there is no special treatment of sex chromosomes, at this point.
}

\seealso{
	\itemize{
	\item \code{\link[qtl]{mqmscan}} - Function called to do single trait analysis
	\item \code{\link[qtl]{scanall}} - Function called to analytse multiple traits in the cross object
	}
}

\examples{
qtl <- c(3,15,3,7)							# QTL at chromosome 3
data(map10)								# Mouse genome
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01) # Simulate a Cross (just autosomes)
crossaug <- mqmaugment(cross,neglect=1)						# Augmentation of the data
result <- mqmscan(crossaug)						# Do a MQM scan of the genome
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
