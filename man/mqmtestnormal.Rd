\name{mqmtestnormal}
\alias{mqmtestnormal}

\title{ Pearsons normality test used for MQM }

\description{
  Wraps a pearsons normality test from the nortest package. 
  This function is used in MQM to test the normality of the
  trait under investigation
}

\usage{
mqmtestnormal(cross, pheno.col = 1)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  \item{pheno.col}{ 
Column number in the phenotype matrix which should be used as the phenotype. This can be a vector of integers. % \phenocol
  }
}

\value{
	Boolean indicating normality of the trait in pheno.col
}

\author{ 
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
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

\keyword{ utilities }
