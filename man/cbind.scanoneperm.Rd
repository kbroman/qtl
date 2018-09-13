\name{cbind.scanoneperm}
\alias{cbind.scanoneperm}

\title{Combine columns from multiple scanone permutation results}

\description{
  Concatenate the columns from different runs of
  \code{\link{scanone}} with \code{n.perm > 0}.
}

\usage{
\method{cbind}{scanoneperm}(\dots, labels)
}

\arguments{
 \item{\dots}{A set of objects of class \code{scanoneperm}. These are
   the permutation results from \code{\link{scanone}} (that is,
   when \code{n.perm > 0}), generally run with different phenotypes or
   methods.}
 \item{labels}{A vector of character strings, of length 1 or of the same
   length as the input \code{\dots}, to be appended to the column names
   in the output.}
}

\value{The concatenated input, as a \code{scanoneperm} object. If
  different numbers of permutation replicates were used, those columns
  with fewer replicates are padded with missing values (\code{NA}).}

\details{
  The aim of this function is to concatenate the results from multiple
  runs of a permutation test \code{\link{scanone}}, generally for
  different phenotypes and/or methods, to be used in parallel with
  \code{\link{c.scanone}}. 
}

\author{Karl W Broman, \email{broman@wisc.edu} }

\examples{
data(fake.f2)
fake.f2 <- calc.genoprob(fake.f2)
\dontshow{fake.f2 <- subset(fake.f2, chr=c(18:19,"X"))}
operm1 <- scanone(fake.f2, method="hk", n.perm=10, perm.Xsp=TRUE)
operm2 <- scanone(fake.f2, method="em", n.perm=5, perm.Xsp=TRUE)

operm <- cbind(operm1, operm2, labels=c("hk","em"))
summary(operm)
}

\seealso{ \code{\link{summary.scanoneperm}},
  \code{\link{scanone}}, \code{\link{c.scanoneperm}},
  \code{\link{c.scanone}} } 

\keyword{manip}
