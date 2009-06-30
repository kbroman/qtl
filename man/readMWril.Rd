\name{readMWril}
\alias{readMWril}

\title{Read data for 4- or 8-way RIL}

\description{
  Data for a set of 4- or 8-way recombinant inbred lines (RIL) is read
  from a pair of comma-delimited files and converted 
  into an object of class \code{cross}.
  We require chromosome
  assignments for the genetic markers, and assume that markers are in
  their correct order.
}

\usage{
readMWril(dir="", rilfile, founderfile,
          type=c("ri4self", "ri4sib", "ri8self", "ri8sib"),
          na.strings=c("-","NA"), rotate=FALSE, \dots)
}

\arguments{
  \item{dir}{Directory in which the data files will be found.  In
    Windows, use forward slashes (\code{"/"}) or double backslashes
    (\code{"\\\\"}) to specify directory trees.}
  \item{rilfile}{Comma-delimited file for the RIL, in the \code{"csv"}
    format described in the help file for \code{\link[qtl]{read.cross}}.}
  \item{founderfile}{File with founder strains' genotypes, in the same
  orientation as the \code{rilfile}, but with just marker names and the
  founders' marker genotypes.}
  \item{type}{The type of RIL.}
  \item{na.strings}{A vector of strings which are to be interpreted as
    missing values.  For the
    \code{csv} formats, these are interpreted globally
    for the entire
    file, so missing value codes in phenotypes must not be valid
    genotypes, and vice versa. For the \code{gary} format, these are
    used only for the phenotype data.} 
  \item{rotate}{If TRUE, the \code{rilfile} and \code{founderfile} are
  rotated (really transposed), with rows corresponding to markers
  and columns corresponding to individuals.}
  \item{\dots}{Additional arguments, passed to the function
    \code{\link[utils]{read.table}} in the case of 
    \code{csv} and \code{csvr} formats.  In particular, one may use the
    argument
    \code{sep} to specify the field separator (the default is a comma)
    and \code{dec} to specify the character used for the decimal point
    (the default is a period).}
}

\value{
  An object of class \code{cross}; see the help file for
  \code{\link[qtl]{read.cross}} for details.

  An additional component \code{crosses} is included; this is a matrix
  indicating the crosses used to generate the RIL.
}

\details{
  The \code{rilfile} should include a phenotype \code{cross} containing
  character strings of the form \code{ABCDEFGH}, indicating the cross
  used to generate each RIL.  

  The founder strains in the \code{founderfile} should be the strains
  \code{A}, \code{B}, \code{C}, \dots, as indicated in the \code{cross}
  phenotype.  
}



\examples{
\dontrun{
ril <- read.cross("../Data", "ril_data.csv", "founder_geno.csv", "ri4self",
                  rotate=TRUE)}
} 

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link[qtl]{read.cross}}, \code{\link[qtl]{sim.cross}} }

\keyword{IO}
