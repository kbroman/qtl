\name{mqmloadanalysis}
\alias{mqmloadanalysis}

\title{Enables loading of previously saved MQM analysis from file} 

\description{
  Enables loading of previously saved MQM analysis from file
}

\usage{
mqmloadanalysis(file = "MQM_output.txt")
}

\arguments{
  \item{file}{ File to load MQM analysis from }
}

%\details{
%	Not available yet
%}

\value{
Returns a matrix with QTLscores for each location on the chromosome, which can be plotted using the \code{\link{plot.scanone}} function.
}

\author{ 
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
}

\seealso{ 
  \itemize{
    % \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item MQM-tutorial - The online MQM tutorial
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Find the most informative cofactors through backward elimination
    \item \code{\link{mqmsetcofactors}} - Set cofactors
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
\dontrun{
qtl <- c(3,15,1,0)					# QTL at chromosome 3
data(map10)						# Mouse genetic map
cross <- sim.cross(map10,qtl,n=100)			# Simulate a Cross
result <- mqmscan(cross,file="aaa.txt")			# Do a MQM scan of the genome and save it to aaa.txt
res_reloaded <- mqmloadanalysis("aaa.txt")		#Reads in and plots the file aaa.txt
}
}

\keyword{ IO }
