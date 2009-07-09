\name{MQMCofactorsEach}
\alias{MQMCofactorsEach}

\title{Create cofactors to be used with scanMQM} 

\description{
  Return a list of cofactors ready to be used with the scanMQM function
}

\usage{
MQMCofactorsEach(cross, each = 3, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}} for details.} 
  \item{each}{ Every 'each' marker will be used as a cofactor }
  \item{verbose}{If TRUE, print tracing information.}
}

%\details{
%	Not available yet
%}

\value{
	An formatted list of cofactors to be used with the \code{\link[qtl]{scanMQM}} algorithm
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} }

\seealso{ 
\itemize{
	\item \code{\link[qtl]{MQM}} - MQM main information
	\item \code{\link[qtl]{scanMQM}} - Single trait analysis
	\item \code{\link[qtl]{scanall}} - Analyse multiple traits 
}
}

\examples{
qtl <- c(3,15,3,7)						# QTL at chromosome 3
data(map10)							# Mouse genome
\dontshow{map10 <- map10[1:5]}
cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
a <- MQMCofactorsEach(cross,3)					# Set Cofactors each third marker
result <- scanMQM(cross,a)					# Do a backward model selection
}

\keyword{ models }
