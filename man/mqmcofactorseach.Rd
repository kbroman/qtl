\name{mqmsetcofactors}
\alias{mqmsetcofactors}

\title{Create cofactors to be used with mqm} 

\description{
  Return a list of cofactors ready to be used with the mqm function
}

\usage{
mqmsetcofactors(cross, each = 3, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.} 
  \item{each}{ Every 'each' marker will be used as a cofactor }
  \item{verbose}{If TRUE, print tracing information.}
}

\value{
	An formatted list of cofactors to be used with the $cofactors$ parameter of the \code{\link{mqm}} algorithm.
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} }

\seealso{ 
  \itemize{
  	\item \code{\link{MQM}} - MQM main information
  	\item \code{\link{mqm}} - Single / multiple trait QTL analysis using MQM
  }
}

\examples{
qtl <- c(3,15,3,7)						# QTL at chromosome 3
data(map10)							# Mouse genome
\dontshow{map10 <- map10[1:5]}
cross <- sim.cross(map10,qtl,n=100)				# Simulate a Cross
a <- mqmsetcofactors(cross,3)					# Set Cofactors each third marker
result <- mqm(cross,a)					# Do a backward model selection
}

\keyword{ models }
