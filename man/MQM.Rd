\name{mqmintro}
\alias{mqmintro}
\alias{MQM}

\title{Multiple QTL Model (MQM) mapping}

\description{
Overview of the MQM mapping functions
}

\usage{
}

\section{Introduction}{
MQM is an implementation of multiple QTL model mapping as described by
Ritsert Jansen (see references below). A genome scan is performed by
moving a QTL along the chromosomes while using a preselected set of
markers as cofactors. If the QTL gets close to a cofactor the cofactor
gets dropped. The backward elimination procedure starts with a
multiple-regression model on an augmented dataset, using a full set of
cofactors evenly spread over the genome. A maximum likelihood (ML)
approach is used for cofactor selection.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\note{
% \input{"inst/docs/mqm/limitations.txt"}

At the moment, MQM is limited to experimental crosses F2, RIL and BC.
Also MQM does not treat sex chromosomes differently from autosomal
chromosomes - though one can introduce sex as a cofactor. Another
limitation is that MQM handles one phenotype at a time. Finally,
automatic backward elimination currently does not consider the marker
density at locations. Future versions of MQM may improve on this.

% -----^^ inst/docs/mqm/limitations.txt ^^-----
}

\seealso{
	\itemize{
% \input{"inst/docs/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning

% -----^^ inst/docs/mqm/standard_seealso.txt ^^-----
	}
}

\examples{
% \input{"inst/docs/mqm/standard_example.txt"}
data(map10)                    # Mouse dataset
qtl <- c(3,15,3,7)             # QTL location

# simulate a cross (autosomes)
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01)

# MQM
crossaug <- mqmaugment(cross)  # Augmentation
result <- mqmscan(crossaug)    # Scan

# show LOD interval of the QTL on chr 3
lodint(result,chr=3)

cat(result$extra$Nind,'real individuals retained in dataset', 
    result$extra$Naug,'individuals augmented')

% -----^^ inst/docs/mqm/standard_example.txt ^^-----
}

\references{
% \input{"inst/docs/mqm/standard_references.txt"}
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
% -----^^ inst/docs/mqm/standard_references.txt ^^-----
}

\keyword{models}
