\name{c.scanone}
\alias{c.scanone}
\alias{cbind.scanone}

\title{Combine columns from multiple scanone results}

\description{
  Concatenate the columns from different runs of
  \code{\link{scanone}}.
}

\usage{
\method{c}{scanone}(\dots, labels)
\method{cbind}{scanone}(\dots, labels)
}

\arguments{
 \item{\dots}{A set of objects of class \code{scanone}. (This can also
 be a list of \code{scanone} objects.) These are
   the results from  \code{\link{scanone}} (with \code{n.perm=0}),
   generally run with different phenotypes or methods.  All must 
   conform with each other, meaning that
   \code{\link{calc.genoprob}} and/or \code{\link{sim.geno}}
   were run with the same values for \code{step} and \code{off.end} and
   with data having the same genetic map.}
 \item{labels}{A vector of character strings, of length 1 or of the same
   length as the input, to be appended to the column names
   in the output.}
}

\value{The concatenated input, as a \code{scanone} object.}

\details{
  The aim of this function is to concatenate the results from multiple
  runs \code{\link{scanone}}, generally for
  different phenotypes and/or methods, to be used in parallel with
  \code{\link{summary.scanone}}. 
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=18:19)}
out.hk <- scanone(fake.f2, method="hk")
out.np <- scanone(fake.f2, model="np")

out <- c(out.hk, out.np, labels=c("hk","np"))
plot(out, lod=1:2, col=c("blue", "red"))
}

\seealso{ \code{\link{summary.scanone}},
  \code{\link{scanone}}, \code{\link{cbind.scanoneperm}} } 

\keyword{manip}
