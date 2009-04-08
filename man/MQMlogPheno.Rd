\name{MQMlogPheno}
\alias{MQMlogPheno}

\title{ 
	MQMlogPheno - Helper function to logtransform phenotypesin R/QTL cross objects 
}
\description{
  Helper function to logtransform phenotypes in an R/QTL cross objects cross&pheno <- log(cross$pheno)
}
\usage{
MQMlogPheno(cross = NULL, Phenot = NULL)
}

\arguments{
  \item{cross}{ Cross object from R/QTL }
  \item{Phenot}{ Phenot to logtransform (if none is specified all phenotypes are transformed) }
}
\details{
	Not available yet
}
\value{
	Cross object with log transformed phenotypes (specified bij Phenot, or all are transformed when Phenot=NULL) 
}

\author{  Danny Arends \email{Danny.Arends@gmail.com} }

\note{ 
	If u find anything wrong with the algorithm please notify me at: \email{Danny.Arends@gmail.com}
}
\seealso{ 
\itemize{
	\item \code{\link{scanMQM}} - Function called to do single trait analysis
	\item \code{\link{scanMQMall}} - Function called to do multiple trait analysis
}
}
\examples{
	data(multitrait)
	plot(scanMQM(MQMaugment(multitrait)))
	plot(scanMQM(MQMaugment(multitrait),doLOG=1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ QTL }
\keyword{ Mapping }
\keyword{ Selection }
