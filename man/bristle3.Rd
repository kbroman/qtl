\name{bristle3}
\alias{bristle3}
\docType{data}

\title{Data on bristle number in Drosophila}

\description{
  Data from bristle number in chromosome 3 recombinant isogenic lines of
  \emph{Drosophila melanogaster}.
}

\usage{data(bristle3)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details.
}

\details{
  There are 66 chromosome 3 recombinant isogenic lines, derived from
  inbred lines that were selected for low (A) and high (B) abdominal
  bristle numbers.  A recombinant chromosome 3 was placed in an isogenic
  low background.

  There are eight phenotypes: the average and SD of the number of
  abdominal and sternopleural bristles in males and females for each
  line.

  Each line is typed at 29 genetic markers on chromosome 3.
}

\references{
  Long, A. D., Mullaney, S. L., Reid, L. A., Fry, J. D., Langley,
  C. H. and MacKay, T. F. C. (1995) High resolution mapping of genetic
  factors affecting abdominal bristle number in \emph{Drosophila
    melanogaster}.  \emph{Genetics} \bold{139}, 1273--1291.
}

%\source{
%}

\examples{
data(bristle3)
# Summaries
summary(bristle3)
plot(bristle3)

# genome scan for each of the average phenotypes
bristle3 <- calc.genoprob(bristle3, step=2)
out <- scanone(bristle3, pheno.col=c(1,3,5,7))

# Plot the results
    # maximum LOD score among four phenotypes
ym <- max(apply(out[,-(1:2)], 2, max))
plot(out, lod=1:3, ylim=c(0,ym))
plot(out, lod=4, add=TRUE, col="green")
}

\seealso{ \code{\link{bristleX}}, \code{\link{listeria}},
  \code{\link{fake.bc}}, \code{\link{fake.f2}},
  \code{\link{fake.4way}}, \code{\link{hyper}} }

\keyword{datasets}
