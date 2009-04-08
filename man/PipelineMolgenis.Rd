\name{PipelineMolgenis}
\alias{PipelineMolgenis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ PipelineMolgenis - Automatic MQM-QTL analysis for Molgenis database systems }
\description{
	Automatic MQM-QTL analysis for Molgenis database systems
}
\usage{
PipelineMolgenis(DBmarkerID, DBtraitID, name = "MQMResults", DBpath, each = 0, n.clusters = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DBmarkerID}{ Database ID of the matrix with individuals*markers }
  \item{DBtraitID}{ Database ID of the matrix with individuals*traits }
  \item{name}{ Name of the table to store the results in }
  \item{DBpath}{ Path to the molgenis database }
  \item{each}{ Every 'each' marker will be used as a cofactor, for dense markersets please set each > 2. }
  \item{n.clusters}{ Number of child processes to split the job into }
  \item{\dots}{ Parameters passed on to \link{scanMQM}  }
}
\details{
  Mostly provided as an example on how-to use use parrallel computation in building a pipeline.
	The Molgenispipeline is comprised of two main parts:
	\itemize{
		\item 1A: Data retrieval from a molgenis database 
		\item 1B: Formatting into a \code{cross} object. See \code{\link[qtl]{read.cross}} for details. 
		\item 2A.1: MQM QTL Model selection, with specified cofactors
		\item 2A.2: MQM QTL mapping using the selected model
		\item 2B: Storage of QTLmap into the database
	}
  When one has the library SNOW installed the second main part is run in parrallel (see the parameter \code{n.clusters}). Otherwise the traits are 
}
\value{
  No values returned ;) (they are stored in the database)
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; M.A. Swertz \email{m.a.swertz@rug.nl} ; R.C. Janssen \email{r.c.jansen@rug.nl} }

\note{  
  The library RCurl from bioconductor should be installed to use a molgenis interface.
  The library SNOW should be installed to use a parrallel MQM analysis. 
}
\seealso{ 
\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do single trait analysis
	\item \code{\link{ResultsToMolgenis}} - Function called to store results to a molgenis database
	\item \code{\link{PipelineMolgenis}} - Automatic MQM-QTL analysis for Molgenis database systems	
}
}
\examples{
	#Setup variables used
	DBmarkerID=298
	DBtraitID=181
	name="Test"
	DBpath="http://celtic.service.rug.nl:8080/molgenis4rsandbox"
	#Call the pipeline
	PipelineMolgenis(DBmarkerID,DBtraitID,name,DBpath,each=3)
}
\references{
	\url{http://www.molgenis.org} \& \url{http://www.Xgap.org} \& \url{http://www.omegahat.org/RCurl/}
	\itemize{
		\item Temple Lang D. (2008), R as a Web Client – the RCurl package \emp{Journal of Statistical Software}
		\item Swertz MA, Jansen RC. (2007), Beyond standardization: dynamic software infrastructures for systems biology. \emph{Nat Rev Genet.} \bold{3}, 235--243.
		\item Swertz MA, De Brock EO, Van Hijum SA, De Jong A, Buist G, Baerends RJ, Kok J, Kuipers OP, Jansen RC. (2004), Molecular Genetics Information System (MOLGENIS): alternatives in developing local experimental genomics databases. \emph{Bioinformatics},\bold{13}, 2075--2083.
	}
}
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
