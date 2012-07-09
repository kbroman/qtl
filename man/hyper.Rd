\name{hyper}
\alias{hyper}
\docType{data}

\title{Data on hypertension}

\description{ 
  Data from an experiment on hypertension in the mouse.
} 

\usage{data(hyper)}

\format{
  An object of class \code{cross}.  See \code{\link{read.cross}}
  for details. 
}  

\details{
There are 250 male backcross individuals typed at 174 markers (actually one
contains only missing values), including 4 on the X chromosome, with one
phenotype. 

The phenotype is the blood pressure.  See the reference below.  Note
that, for most markers, genotypes are available on only the individuals
with extreme phenotypes.  At many markers, only recombinant individuals
were typed.
}


\references{
  Sugiyama, F., Churchill, G. A., Higgens, D. C., Johns, C.,
  Makaritsis, K. P., Gavras, H. and Paigen, B. (2001) Concordance of
  murine quantitative trait loci for salt-induced hypertension with rat
  and human loci. \emph{Genomics} \bold{71}, 70--77. 
}
  
\source{
  Bev Paigen and Gary Churchill (The Jackson Laboratory, Bar Harbor,
  Maine)
  \url{http://www.jax.org/research/churchill/datasets/qtl/qtlarchive}   
}

\examples{
data(hyper)
summary(hyper)
plot(hyper)

# Note the selective genotyping
\dontrun{plotMissing(hyper, reorder=TRUE)}

# A marker on c14 has no data; remove it
hyper <- drop.nullmarkers(hyper)
}

\seealso{ \code{\link{fake.bc}}, \code{\link{fake.f2}},
  \code{\link{fake.4way}}, \code{\link{listeria}},
  \code{\link{bristle3}}, \code{\link{bristleX}} } 

\keyword{datasets}
