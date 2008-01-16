\name{sim.cross}
\alias{sim.cross}

\title{Simulate a QTL experiment}

\description{
  Simulates data for a QTL experiment using a model in which QTLs act additively.
}

\usage{
sim.cross(map, model=NULL, n.ind=100, type=c("f2", "bc", "4way"),
          error.prob=0, missing.prob=0, partial.missing.prob=0,
          keep.qtlgeno=TRUE, keep.errorind=TRUE, m=0, p=0,
	  map.function=c("haldane","kosambi","c-f","morgan"))
}

\arguments{
 \item{map}{A list whose components are vectors containing the marker
   locations on each of the chromosomes.}
 \item{model}{A matrix where each row corresponds to a
 different QTL, and gives the chromosome number, cM position and
 effects of the QTL.}
 \item{n.ind}{Number of individuals to simulate.}
 \item{type}{Indicates whether to simulate an intercross (\code{f2}),
   a backcross (\code{bc}), or a phase-known 4-way cross (\code{4way}).}
 \item{error.prob}{The genotyping error rate.}
 \item{missing.prob}{The rate of missing genotypes.}
 \item{partial.missing.prob}{When simulating an intercross or 4-way
   cross, this gives the rate at which markers will be incompletely
   informative (i.e., dominant or recessive).}
 \item{keep.qtlgeno}{If TRUE, genotypes for the simulated QTLs will be
   included in the output.}
 \item{keep.errorind}{If TRUE, and if \code{error.prob > 0}, the
   identity of genotyping errors will be included in the output.}
 \item{m}{Interference parameter; a non-negative integer.  0 corresponds
   to no interference.}
 \item{p}{Probability that a chiasma comes from the no-interference
   mechanism} 
 \item{map.function}{Indicates whether to use the Haldane, Kosambi, 
   Carter-Falconer, or Morgan map function when converting genetic
   distances into recombination fractions.}
}

\details{
  Meiosis is assumed to follow the Stahl model for crossover
  interference (see the references, below), of which the no interference
  model and the chi-square model are special cases.  Chiasmata on the
  four-strand bundle are a superposition of chiasmata from two different
  mechanisms.  With probability \code{p}, they arise by a mechanism
  exhibiting no interference; the remainder come from a chi-square model
  with inteference parameter \code{m}.  Note that \code{m=0} corresponds
  to no interference, and with \code{p=0}, one gets a pure chi-square
  model.
  
  If a chromosomes has class \code{X}, it is assumed to be the X
  chromosome, and is assumed to be segregating in the cross.  Thus, in
  an intercross, it is segregating like a backcross chromosome.  In a
  4-way cross, a second phenotype, \code{sex}, will be generated.

  QTLs are assumed to act additively, and the residual phenotypic
  variation is assumed to be normally distributed with variance 1.

  For a backcross, the effect of a QTL is a single number corresponding
  to the difference between the homozygote and the heterozygote.

  For an intercross, the effect of a QTL is a pair of numbers,
  (\eqn{a,d}), where \eqn{a} is the additive effect (half the difference
  between the homozygotes) and \eqn{d} is the dominance deviation (the
  difference between the heterozygote and the midpoint between the
  homozygotes).

  For a four-way cross, the effect of a QTL is a set of three numbers,
  (\eqn{a,b,c}), where, in the case of one QTL, the mean phenotype,
  conditional on the QTL genotyping being AC, BC, AD or BD, is  \eqn{a},
  \eqn{b}, \eqn{c} or 0, respectively.
}  

\value{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}} for
  details. 

  If \code{keep.qtlgeno} is TRUE, the cross object will contain a
  component \code{qtlgeno} which is a matrix containing the QTL
  genotypes (with complete data and no errors), coded as in the genotype
  data.

  If \code{keep.errorind} is TRUE and errors were simulated, each
  component of \code{geno} will each contain a matrix \code{errors},
  with 1's indicating simulated genotyping errors.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} }

\seealso{ \code{\link[qtl]{sim.map}}, \code{\link[qtl]{read.cross}},
    \code{\link[qtl]{fake.f2}}, \code{\link[qtl]{fake.bc}}
    \code{\link[qtl]{fake.4way}} }

\examples{
# simulate a genetic map
map <- sim.map()

# simulate 250 intercross individuals with 2 QTLs
fake <- sim.cross(map, type="f2", n.ind=250,
                  model = rbind(c(1,45,1,1),c(5,20,0.5,-0.5)))
}
                    
\references{
  Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
  interference in arabidopsis.  \emph{Genetics} \bold{160}, 1631--1639.

  Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
  interference as a function of genetic distance. \emph{Genetics}
  \bold{133}, 681--691.

  Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis
  of crossover interference using the chi-square model.  \emph{Genetics}
  \bold{139}, 1045--1056.

}
  
\keyword{datagen}
