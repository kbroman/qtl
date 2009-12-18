\name{mqmpermutationresult}
\alias{mqmpermutationresult}

\title{Convert mqmmulti objects into a scanoneperm object}

\description{
	Function to convert \code{mqmmulti} objects into a \code{scanoneperm} 
  object, this allows the use of R/qtl methods for permutation analysis
  that do not support the output of a multiple QTL scan using mqm's 
  outputstructure.
}

\usage{
mqmpermutationresult(mqmpermutationresult = NULL)
}

\arguments{
  \item{mqmpermutationresult}{ \code{mqmmulti} object obtained 
    after performing permutations on a single trait.using the 
    function \code{\link[qtl]{mqmpermute}} 
  }
}

%\details{
%
%}

\value{
	Output of the algorithm is a \code{scanoneperm} object. See also: \code{\link[qtl]{summary.scanoneperm}}
}

\author{ Danny Arends \email{Danny.Arends@gmail.com} ; R.C. Jansen \email{r.c.jansen@rug.nl} }

\seealso{ 
  \itemize{
    \item \code{\link[qtl]{mqmintro}} - Introduction into MQM
    \item \code{\link[qtl]{mqmoverview}} - Overview of the MQM specific additions to R/qtl
  	\item \code{\link[qtl]{mqmpermute}} - Estimate confidence interfalls for QTL location using MQM
  	\item \code{\link[qtl]{summary.scanoneperm}} - LOD thresholds from scanone permutation results
  }
}

\examples{
f2qtl <- c(3,15,3,7)					# QTL at chromosome 3
data(map10)						# Mouse genome
\dontshow{map10 <- map10[1:3]}
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")	# Simulate a F2 Cross
\dontrun{f2result <- bootstrap(f2cross,mapfunction=mqm)			# Bootstrap MQM mapping on the f2cross
}\dontshow{f2result <- bootstrap(f2cross,mapfunction=mqm, n.run=2, b.size=2)}
f2perm <- mqmpermutationresult(f2result)			# Create a permutation object
summary(f2perm)					# What LOD score is considered significant ?
}

\keyword{ models }
