\name{ResultsToMolgenis}
\alias{ResultsToMolgenis}

\title{ ResultsToMolgenis - Helperfunction to store results from a MQM-qtl run into a molgenis database }
\description{
  Helperfunction to store results from a MQM-qtl run into a a molgenis database
}
\usage{
ResultsToMolgenis(intervalQTLmap = NULL, name = NULL, DBpath = NULL, Fupdate = 0)
}
\arguments{
  \item{intervalQTLmap}{ Results from either scanone, scanMQM or scanMQMall }
  \item{name}{ Name of the table to store the results in }
  \item{Trait_num}{ Analysis is stored in this column in the database. This parameter is used for \link{PipelineMolgenis} }
  \item{DBpath}{ Path to the molgenis database }
  \item{verbose}{ Print to screen (DEFAULT: TRUE) }
}
\details{
 The library RCurl from bioconductor should be installed to use a molgenis interface.
}
\value{
	No return value
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; M.A. Swertz \email{m.a.swertz@rug.nl} ; R.C. Janssen \email{r.c.jansen@rug.nl} }

\note{
The library RCurl should be installed before attempting to use this function. This library can be obtained from bioconductor or omegahat.
\itemize{
	\item \url{http://www.bioconductor.org} or \url{http://www.omegahat.org/RCurl/}
	\item \url{http://www.molgenis.org}
	\item \url{http://www.Xgap.org}
}
If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}

\seealso{ 
\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{CrossFromMolgenis}} - Function called to store results to a molgenis database
	\item \code{\link{PipelineMolgenis}} - Automatic MQM-QTL analysis for Molgenis database systems		
}
}
\examples{
cMol <- CrossFromMolgenis()			#Make a crossobject from the database
result <- scanMQMall(cMol)			#Basic execution of scanMQMall
ResultsToMolgenis(result,"QTL-mqm")	#Store the result from MQM into molgenis database
}
\references{
	\url{http://www.molgenis.org} \& \url{http://www.Xgap.org} \& \url{http://www.omegahat.org/RCurl/}	
	\itemize{
		\item Temple Lang D. (2008), R as a Web Client – the RCurl package \emp{Journal of Statistical Software}
		\item Swertz MA, Jansen RC. (2007), Beyond standardization: dynamic software infrastructures for systems biology. \emph{Nat Rev Genet.} \bold{3}, 235--243.
		\item Swertz MA, De Brock EO, Van Hijum SA, De Jong A, Buist G, Baerends RJ, Kok J, Kuipers OP, Jansen RC. (2004), Molecular Genetics Information System (MOLGENIS): alternatives in developing local experimental genomics databases. \emph{Bioinformatics},\bold{13}, 2075--2083.
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
