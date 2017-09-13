\name{mqmprocesspermutation}
\alias{mqmprocesspermutation}

\title{Convert mqmmulti objects into a scanoneperm object}

\description{
    Function to convert \code{mqmmulti} objects into a \code{scanoneperm}
  object, this allows the use of R/qtl methods for permutation analysis
  that do not support the output of a multiple QTL scan using mqm's
  outputstructure.
}

\usage{
mqmprocesspermutation(mqmpermutationresult = NULL)
}

\arguments{
  \item{mqmpermutationresult}{ \code{mqmmulti} object obtained
    after performing permutations on a single trait.using the
    function \code{\link{mqmpermutation}}
  }
}

\value{
    Output of the algorithm is a \code{scanoneperm} object. See also: \code{\link{summary.scanoneperm}}
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\seealso{
  \itemize{
% \input{"inst/doc/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{https://rqtl.org/tutorials/MQM-tour.pdf}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Set cofactors using marker density
    \item \code{\link{mqmsetcofactors}} - Set cofactors at fixed locations
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/doc/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
# QTL at chromosome 3
f2qtl <- c(3,15,1,0)

# Mouse genetic map
data(map10)

\dontshow{map10 <- map10[1:3]}
# Simulate a F2 Cross
f2cross <- sim.cross(map10,f2qtl,n=100,type="f2")
\dontrun{# Bootstrap MQM mapping on the f2cross
f2result <- mqmpermutation(f2cross,scanfunction=mqmscan)
}\dontshow{f2result <- mqmpermutation(f2cross,scanfunction=mqmscan, n.perm=2, batchsize=2)}

# Create a permutation object
f2perm <- mqmprocesspermutation(f2result)

# What LOD score is considered significant?
summary(f2perm)
}

\keyword{ models }
