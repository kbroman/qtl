\name{mqmintro}
\alias{mqmintro}
\alias{MQM}
\alias{mqm}

\title{Multiple QTL Model (MQM) mapping}

\description{
Overview of the MQM mapping functions
}

\usage{
}

\section{Introduction}{
% \input{"inst/docs/Sources/MQM/mqm/description.txt"}
Multiple QTL model (MQM) mapping is a powerful QTL mapping method
developed by Ritsert C. Jansen\cite{jansen94}.  MQM is an automatic
two-stage procedure in which, in the first stage, important markers
are selected in multiple regression on markers. In the second stage a
QTL is moved along the chromosomes by using the pre-selected markers
as cofactors, except for the markers flanking the interval under
study. A refined and automated procedure for cases with large numbers
of marker cofactors is included.

MQM for R/qtl is based on the original implementation of MQM mapping
and consists of a three step procedure: (1) data augmentation, (2)
automatic backward model elimination using genetic markers evenly
spread over the genome as cofactors and (3) QTL (interval) mapping
using the most `informative' model through maximum likelihood. The
method internally controls false discovery rates (FDR) and lets users
test different QTL models by elimation of non-significant cofactors. 

% -----^^ inst/docs/Sources/MQM/mqm/description.txt ^^-----
% \input{"inst/docs/Sources/MQM/mqm/advantages_Rd.txt"}
R/qtl-MQM has the following advantages:

\itemize{
  \item Higher power, as long as the QTL explain a reasonable amount of variation
  \item Protection against overfitting, because it uses the residual variance from the full model. For this reason more parameters (cofactors) can be used compared to, for example, CIM
  \item Prevention of ghost QTL (between two QTL in coupling phase)
  \item Detection of negating QTL (QTL in repulsion phase)}
% -----^^ inst/docs/Sources/MQM/mqm/advantages_Rd.txt ^^-----
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\note{
% \input{"inst/docs/Sources/MQM/mqm/limitations.txt"}

The current implementation of R/qtl-MQM has the following recognized
limitations: (1) MQM is limited to experimental crosses F2, RIL and
BC.  (2) MQM augmentation drops entire individuals from the equation
when a lot of marker data is missing.  (3) MQM does not treat sex
chromosomes differently from autosomal chromosomes - though one can
introduce sex as a cofactor. Finally, (4) automatic backward
elimination sets cofactors at fixed locations and ignores underlying
known marker density. Future versions of R/qtl-MQM may improve on
these points. Check the changelog for updates.

% -----^^ inst/docs/Sources/MQM/mqm/limitations.txt ^^-----
}

\seealso{
	\itemize{
    \item The MQM tutorial (online)
% \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
	}
}

\examples{
% \input{"inst/docs/Sources/MQM/mqm/standard_example.txt"}
data(map10)                    # Mouse dataset

# simulate a cross (autosomes)
qtl <- c(3,15,3,7)             # QTL location
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01)

# MQM
crossaug <- mqmaugment(cross)  # Augmentation
cat(crossaug$mqm$Nind,'real individuals retained in dataset',
    crossaug$mqm$Naug,'individuals augmented')

result <- mqmscan(crossaug)    # Scan

# show LOD interval of the QTL on chr 3
lodint(result,chr=3)
% -----^^ inst/docs/Sources/MQM/mqm/standard_example.txt ^^-----
}

\references{
\itemize{
% \input{"inst/docs/Sources/MQM/mqm/standard_references.txt"}
		\item Jansen RC, (2007) Quantitative trait loci in inbred lines. Chapter 18 of \emph{Handbook of Stat. Genetics} 3rd edition. John Wiley & Sons, Ltd.
		\item Jansen RC, Nap JP (2001), Genetical genomics: the added value from segregation. \emph{Trends in Genetics}, \bold{17}, 388--391.
		\item Jansen RC, Stam P (1994), High resolution of quantitative traits into multiple loci via interval mapping. \emph{Genetics}, \bold{136}, 1447--1455.
		\item Jansen RC (1993), Interval mapping of multiple quantitative trait loci. \emph{Genetics}, \bold{135}, 205--211.
		\item Swertz MA, Jansen RC. (2007), Beyond standardization: dynamic software infrastructures for systems biology. \emph{Nat Rev Genet.} \bold{3}, 235--243.
		\item Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977) Maximum likelihood from incomplete data via the EM algorithm.  \emph{J. Roy. Statist. Soc.} B, \bold{39}, 1--38.

% -----^^ inst/docs/Sources/MQM/mqm/standard_references.txt ^^-----
}
}

\keyword{models}
