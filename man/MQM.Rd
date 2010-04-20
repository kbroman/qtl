\name{MQM}
\alias{MQM}

\title{Introduction to Multiple QTL Model (MQM) mapping}

\description{
Overview of the MQM mapping functions
}

%\usage{
%}

\section{Introduction}{
% \input{"inst/docs/Sources/MQM/mqm/description.txt"}
Multiple QTL Model (MQM) mapping provides a sensitive approach for
mapping quantititive trait loci (QTL) in experimental populations. MQM
adds higher statistical power compared to many other methods.  The
theoretical framework of MQM mapping was introduced and explored by
Ritsert Jansen, explained in the `Handbook of Statistical Genetics'
(see references), and used effectively in practical research, with the
commercial `mapqtl' software package. Here we present the first free
and open source implementation of MQM, with extra features like high
performance parallelization on multi-CPU computers, new plots and
significance testing.

MQM is an automatic three-stage procedure in which, in the first
stage, missing data is `augmented'. In other words, rather than guessing one
likely genotype, multiple
genotypes are modeled with their estimated
probabilities.  In the second stage important markers are selected by
multiple regression and backward elimination. In the third stage a QTL is moved
along the chromosomes using these pre-selected markers as cofactors,
except for the markers flanking the interval under study. QTL are
(interval) mapped using the most `informative' model through maximum
likelihood.  A refined and automated procedure for cases with large
numbers of marker cofactors is included.  The method internally
controls false discovery rates (FDR) and lets users test different QTL
models by elimation of non-significant cofactors. 

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

The current implementation of R/qtl-MQM has the following 
limitations: (1) MQM is limited to experimental crosses F2, RIL and
BC, (2) MQM augmentation drops entire individuals from the equation
when a lot of marker data is missing, (3) MQM does not treat sex
chromosomes differently from autosomal chromosomes - though one can
introduce sex as a cofactor. Finally, (4) automatic backward
elimination sets cofactors at fixed locations and ignores underlying
known marker density. Future versions of R/qtl-MQM may improve on
these points. Check the website and changelog for updates.

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
    \item \code{\link{mqmsetcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmpermutation}} - Estimate significance levels
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
