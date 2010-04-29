\name{mqmplot.permutations}
\alias{mqmplot.permutations}

\title{Plot results from mqmpermutation}

\description{
	Plotting routine to display the results from a permutation QTL scan. (the
  output of \code{\link{mqmpermutation}})
}

\usage{
mqmplot.permutations(permutationresult, \dots)
}

\arguments{
  \item{permutationresult}{ \code{mqmmulti} object returned by
  \code{\link{mqmpermutation}} permutation analysis.
  }
  \item{\dots}{ Extra arguments passed to \code{\link{polyplot}}}
  }

%\details{
%  ~~ If necessary, more details than the description above ~~
%}

\value{ No value returned (plotting routine) }

\author{ 
Danny Arends \email{danny.arends@gmail.com} % \dannyauthor
  , Rutger Brouwer
}

\seealso{ 
\itemize{
  % \input{"inst/docs/Sources/MQM/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item MQM-tutorial - The online MQM tutorial
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Set cofactors using marker density
    \item \code{\link{mqmsetcofactors}} - Set cofactors at fixed locations
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
#Simulated F2 Population
f2qtl <- c(3,15,1,0)                                    # QTL at chromosome 3
data(map10)                                             # Mouse genetic map
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")       # Simulate a F2 Cross
f2cross <- calc.genoprob(f2cross)
\dontrun{f2result <- mqmpermutation(f2cross, n.perm=1000, method="permutation")	# Permutations to obtain significance threshold
}\dontshow{f2result <- mqmpermutation(f2cross, n.perm=2, batchsize=2, method="permutation")}
mqmplot.permutations(f2result)					# Plot results
}
\keyword{ hplot }
