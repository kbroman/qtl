\name{sim.cross}
\alias{sim.cross}

\title{Simulate a QTL experiment}

\description{
  Simulates data for a QTL experiment using a model in which QTLs act additively.
}

\usage{
sim.cross(map, model=NULL, n.ind=100, 
         type=c("f2", "bc", "4way", "risib", "riself",
           "ri4sib", "ri4self", "ri8sib", "ri8self", "bcsft"),
          error.prob=0, missing.prob=0, partial.missing.prob=0,
          keep.qtlgeno=TRUE, keep.errorind=TRUE, m=0, p=0,
	  map.function=c("haldane","kosambi","c-f","morgan"),
          founderGeno, random.cross=TRUE, \dots)
}

\arguments{
 \item{map}{A list whose components are vectors containing the marker
   locations on each of the chromosomes.}
 \item{model}{A matrix where each row corresponds to a
 different QTL, and gives the chromosome number, cM position and
 effects of the QTL.}
 \item{n.ind}{Number of individuals to simulate.}
 \item{type}{Indicates whether to simulate an intercross (\code{f2}),
   a backcross (\code{bc}), a phase-known 4-way cross (\code{4way}),
   or recombinant inbred lines (by selfing or by sib-mating, and with
   the usual 2 founder strains or with 4 or 8 founder strains).}
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
 \item{founderGeno}{For 4- or 8-way RIL, the genotype data of the
   founder strains, as a list whose components are numeric matrices
   (no. markers x no. founders), one for each chromosome.}
 \item{random.cross}{For 4- or 8-way RIL, indicates whether the order of the
    founder strains should be randomized, independently for each RIL, or
    whether all RIL be derived from a common cross.  In the latter case,
    for a 4-way RIL, the cross would be (AxB)x(CxD).}
 \item{\dots}{For \code{type = "bcsft"}, additional arguments passed to \code{sim.cross.bcsft}.}
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

\section{Recombinant inbred lines}{
  In the simulation of recombinant inbred lines (RIL), we simulate a
  single individual from line, and no phenotypes are simulated (so the
  argument \code{model} is ignored).

  The types \code{riself} and \code{risib} are the usual two-way RIL.

  The types \code{ri4self}, \code{ri4sib}, \code{ri8self}, and
  \code{ri8sib} are RIL by selfing or sib-mating derived from four or
  eight founding parental strains.

  For the 4- and 8-way RIL, one must include the genotypes of the
  founding individuals; these may be simulated with
  \code{\link{simFounderSnps}}.  Also, the output cross will
  contain a component \code{cross}, which is a matrix with rows
  corresponding to RIL and columns corresponding to the founders, 
  indicating order of the founder strains in the crosses used to
  generate the RIL.

  The coding of genotypes in 4- and 8-way RIL is rather complicated.  It
  is a binary encoding of which founder strains' genotypes match the
  RIL's genotype at a marker, and not that this is specific to the order
  of the founders in the crosses used to generate the RIL.  For example,
  if an RIL generated from 4 founders has the 1 allele at a SNP, and the
  four founders have SNP alleles 0, 1, 0, 1, then the RIL allele matches
  that of founders B and D.  If the RIL was derived by the cross (AxB)x(CxD),
  then the RIL genotype would be encoded \eqn{2^{2-1} + 2^{3-1} = 6}{2^(2-1) + 2^(3-1) = 6}.
  If the cross was derived by the cross (DxA)x(CxB), then the RIL
  genotype would be encoded \eqn{2^{1-1} + 2^{4-1} = 9}{2^(1-1) + 2^(4-1) = 6}.
  These get reorganized after calls to \code{\link{calc.genoprob}},
  \code{\link{sim.geno}}, or \code{\link{argmax.geno}}, and
  this approach simplifies the hidden Markov model (HMM) code.

  For the 4- and 8-way RIL, genotyping errors are simulated only if the
  founder genotypes are 0/1 SNPs.
}

\value{
  An object of class \code{cross}.  See \code{\link{read.cross}} for
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

\seealso{ \code{\link{sim.map}}, \code{\link{read.cross}},
    \code{\link{fake.f2}}, \code{\link{fake.bc}}
    \code{\link{fake.4way}}, \code{\link{simFounderSnps}} }

\examples{
# simulate a genetic map
map <- sim.map()


### simulate 250 intercross individuals with 2 QTLs
fake <- sim.cross(map, type="f2", n.ind=250,
                  model = rbind(c(1,45,1,1),c(5,20,0.5,-0.5)))


### simulate 100 backcross individuals with 3 QTL
# a 10-cM map model after the mouse
data(map10)

fakebc <- sim.cross(map10, type="bc", n.ind=100,
                    model=rbind(c(1,45,1), c(5,20,1), c(5,50,1)))


### simulate 8-way RIL by sibling mating
# get lengths from the above 10-cM map
L <- ceiling(sapply(map10, max))

# simulate a 1 cM map
themap <- sim.map(L, n.mar=L+1, eq.spacing=TRUE)

# simulate founder genotypes
pg <- simFounderSnps(themap, "8")

# simulate the 8-way RIL by sib mating (256 lines)
ril <- sim.cross(themap, n.ind=256, type="ri8sib", founderGeno=pg)
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

  Broman, K. W. (2005) The genomes of recombinant inbred lines
  \emph{Genetics} \bold{169}, 1133--1146.

  Teuscher, F. and Broman, K. W. (2007) Haplotype probabilities for
  multiple-strain recombinant inbred lines.  \emph{Genetics} \bold{175},
  1267--1274. 
}
  
\keyword{datagen}
