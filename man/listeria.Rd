\name{listeria}
\alias{listeria}
\docType{data}

\title{Data on Listeria monocytogenes susceptibility}

\description{
  Data from an experiment on susceptibility to \emph{Listeria
  monocytogenes} infection in the mouse.
}

\usage{data(listeria)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details.
}

\details{
  There are 120 F2 individuals typed at 133 markers, including 2 on the
  X chromosome, with one phenotype.

  The phenotype is the survival time (in hours) following infection.
  Mice with phenotype 264 hours may be considered to have recovered
  from the infection. See the references below.
}

\references{
  Boyartchuk, V. L., Broman, K. W., Mosher, R. E., D'Orazio
  S. E. F., Starnbach, M. N. and Dietrich, W. F. (2001) Multigenic
  control of \emph{Listeria monocytogenes} susceptibility in
  mice. \emph{Nature Genetics} \bold{27}, 259--260.

  Broman,  K. W. (2003) Mapping quantitative trait loci in the case of a
  spike in the phenotype distribution. \emph{Genetics} \bold{163},
  1169--1175.
}

\source{
  Victor Boyartchuk and William Dietrich (Department of Genetics,
  Harvard Medical School and Howard Hughes Medical Institute)
}

\examples{
data(listeria)
\dontshow{listeria <- subset(listeria,chr=c(1,5,6,13,15), ind=!is.na(listeria$pheno[,1]))}
# Summaries
summary(listeria)
plot(listeria)

# Take log of phenotype
listeria$pheno[,1] <- log2(listeria$pheno[,1])
plot(listeria)

# Genome scan with a two-part model, using log survival
listeria <- calc.genoprob(listeria, step=2)
out <- scanone(listeria, model="2part", method="em",
               upper=TRUE)

# Summary of the results
summary(out, thr=c(5,3,3), format="allpeaks")

# Plot LOD curves for interesting chromosomes
#     (The two-part model gives three LOD scores)
plot(out, chr=c(1,5,6,13,15), lodcolumn=1:3,
     lty=1, col=c("black","red","blue"))
}

\seealso{ \code{\link{fake.bc}}, \code{\link{fake.f2}},
  \code{\link{fake.4way}}, \code{\link{hyper}},
  \code{\link{bristle3}}, \code{\link{bristleX}} }

\keyword{datasets}
