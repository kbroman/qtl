\name{MQMpermObject}
\alias{MQMpermObject}

\title{ MQMpermObject - Function to convert MQMmulti objects into a scanoneperm object from R/QTL }
\description{
	Function to convert MQMmulti objects into a scanoneperm object from R/QTL, this allows us to use the R/QTL methods of permutation analysis
}
\usage{
MQMpermObject(MQMbootresult = NULL)
}

\arguments{
  \item{MQMbootresult}{ MQMmulti object obtained by bootstrapMQM }
}
\details{

}
\value{
	Output of the algorithm is a scanoneperm object See also: \link[qtl]{summary.scanoneperm}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Janssen \email{r.c.jansen@rug.nl} }
\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link{bootstrapMQM}} - Function called to estimate confidence interfalls for MQM
	\item \link[qtl]{summary.scanoneperm} - Summary of a QTL object
}
}
\examples{
library(MQMpackage)
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
f2result <- bootstrapMQM(f2cross)			# Bootstrap MQM mapping on the f2cross
f2perm <- MQMpermObject(f2result)			# Create a permutation object
summary(f2perm)					# What LOD score is considered significant ?

#EXTRA:
numresults <- dim(f2result[[1]])[1]
lines(x=1:numresults,y=rep(summary(f2perm)[1,1],numresults),col="green",lwd=2,lty=2)	# Green line for A=0.05
lines(x=1:numresults,y=rep(summary(f2perm)[2,1],numresults),col="blue",lwd=2,lty=2)		# Blue line for A=0.10
summary(f2perm)										# What LOD score is considered significant ?
summary(f2result[[1]],alpha=0.05,perms=f2perm)		# Get significant QTL's (The first run is the original phenotype original)
}

\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
