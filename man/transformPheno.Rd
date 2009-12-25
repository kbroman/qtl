\name{transformPheno}
\alias{transformPheno}

\title{Transform phenotypes in a cross object}

\description{
  Transform phenotypes in a cross object; by default use a Log
  transformation, though any function may be used.
}

\usage{
transformPheno(cross, pheno.col=1, transf=log, ...)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}} for details.}
 \item{pheno.col}{A vector of numeric indices or
   character strings (indicating phenotypes by name) of phenotypes to be
   transformed.}
 \item{transf}{The function to use in the transformation.}
 \item{\dots}{Additional arguments, to be passed to \code{transf}.}
}

%\details{
%  The selected phenotype columns are transformed via the selected
%  transformation.  
%}

\value{
  The input cross object with the transformed phenotypes
}

\author{  Danny Arends, \email{Danny.Arends@gmail.com}, Karl W. Broman,
  \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{transf}}, \code{\link{mqmscan}}, \code{\link{scanone}} }

\examples{
data(multitrait)

# Log transform all phenotypes
multitrait.log <- transformPheno(multitrait, pheno.col=1:nphe(multitrait))

# Square-root transform all phenotypes
multitrait.sqrt <- transformPheno(multitrait, pheno.col=1:nphe(multitrait),
                                  transf=sqrt)
}

\keyword{ utilities }
