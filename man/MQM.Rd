\name{MQM}
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
    \item \code{\link[qtl]{MQM}} - MQM description and references
    \item \code{\link[qtl]{mqmscan}} - Main MQM single trait analysis
    \item \code{\link[qtl]{mqmall}} - Parallellized traits analysis 
    \item \code{\link[qtl]{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link[qtl]{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link[qtl]{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link[qtl]{scanone}} - Single QTL scanning

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

cat(result$extra$Nind,'real individuals retained in dataset', \
    result$extra$Naug,'individuals augmented')

% -----^^ inst/docs/mqm/standard_example.txt ^^-----
}

\references{
% \input{"inst/docs/mqm/standard_references.txt"}
}

\keyword{models}
