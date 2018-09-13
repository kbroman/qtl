\name{convert.scantwo}
\alias{convert.scantwo}

\title{Convert output from scantwo for R/qtl version 1.03 and earlier}

\description{
  Convert the output from scantwo from the format used in R/qtl version
  1.03 and earlier to that used in version 1.04 and later.
}

\usage{
\method{convert}{scantwo}(object, \dots)
}

\arguments{
 \item{object}{Output from the function \code{\link{scantwo}}, for
   R/qtl version 1.03 and earlier.}
 \item{\dots}{Ignored at this point.}
}

\value{
  The same scanone output, but revised for use with R/qtl version 1.03
  and later.
}

\details{
  Previously, the output from \code{\link{scantwo}} contained the
  full and interaction LOD scores.  In R/qtl version 1.04 and later,
  the output contains the LOD scores from the full and
  additive QTL models.
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
\dontrun{out2.new <- convert(out2.old)}
}

\seealso{ \code{\link{scantwo}}, \code{\link{convert.scanone}} }

\keyword{manip}
