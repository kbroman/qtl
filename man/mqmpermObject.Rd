\name{mqmpermObject}
\alias{mqmpermObject}

\title{Function to convert mqmmulti objects into a scanoneperm object}

\description{
	Function to convert mqmmulti objects into a \code{scanoneperm} object, this allows the use of R/QTL methods for permutation analysis that do not support MQM's internal structure.
}

\usage{
mqmpermObject(mqmbootresult = NULL)
}

\arguments{
  \item{mqmbootresult}{ \code{mqmmulti} object obtained by using the function \code{\link[qtl]{bootstrap}} }
}

%\details{
%
%}

\value{
	Output of the algorithm is a \code{scanoneperm} object See also: \code{\link[qtl]{summary.scanoneperm}}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} }

\seealso{ 
\itemize{
	\item \code{\link[qtl]{bootstrap}} - Function called to estimate confidence interfalls for scanone,cim or mqm
	\item \code{\link[qtl]{summary.scanoneperm}} - Summary of a QTL object
}
}

\examples{
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
\dontrun{f2result <- bootstrap(f2cross,Funktie=mqm)			# Bootstrap MQM mapping on the f2cross
}\dontshow{f2result <- bootstrap(f2cross,Funktie=mqm, n.run=2, b.size=2)}
f2perm <- mqmpermObject(f2result)			# Create a permutation object
summary(f2perm)					# What LOD score is considered significant ?
}

\keyword{ models }
