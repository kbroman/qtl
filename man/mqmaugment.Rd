\name{mqmaugment}
\alias{mqmaugment}

\title{MQM augmentation}

\description{
Fills in missing genotypes for multiple QTL model (MQM) mapping. For each
missing or incomplete marker it fills in (or 'augments') all possible
genotypes, thus creating new candidate 'individuals' and calculates their
probability using information on neighbouring markers and recombination
frequencies. When a genotype of an augmented genotype is less likely than the
\code{minprob} parameter it is dropped from the dataset. The \emph{augmented} list
of individuals is returned in a new \code{crossobject}.  
}

\usage{
mqmaugment(cross, maxaugind=60, minprob=0.1, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link{read.cross}}.} 
  \item{maxaugind}{ 
    Maximum number of augmentations per individual.  The default of 60 allows
    for 20 missing markers for an individual in a BC cross and 30 missing
    markers in an F2. When a large number of markers are missing this default
    number is quickly reached.
  }
  \item{minprob}{ 
    Return individuals with augmented genotypes that are above a probability of
    occurring. minprob is a value between 0 and 1. For example a value of 0.5
    will drop all genotypes that are half as likely than the most likely
    genotype (candidate of the individual). The default value of 0.1 will drop
    all genotypes that are less likely of ocurring than 1 in 10, compared
    against the most likely genotype. Use a value of 1 to return a single
    filled in genotype for each individual.
 }
  \item{verbose}{ Verbose output }
}

\value{
  Returns the cross object with augmented individuals (many individuals 
  from the data set will be repeated multiple times). Some individuals 
  may have been dropped completely when the probability falls below 
  \code{minprob}. An added component to the cross object named \code{mqm} contains
  information on exactly which individuals are retained and repeated.  
}  

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\note{
The sex chromosome 'X' is treated like autosomes during augmentation.
With an F2 the sex chromosome is not considered. This will change in
a future version of MQM.

When augmentation runs out of memory, a warning is given and the returned
object should not be used.

With a large amount of missing data, e.g. larger than 5\%, a low \code{minprob}
may cause the routine to drop a lot of indivuals from the dataset and distort
results. Run with \code{verbose=TRUE} to verify how many individuals are
retained in the dataset or check the return \code{cross$mqm} values. Compare
results by using \code{minprob=1}.
}

\seealso{
  \itemize{
	\item \code{\link{fill.geno}} - Alternative routine for estimating missing data	
% \input{"inst/docs/mqm/standard_seealso.txt"}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmall}} - Parallellized traits analysis 
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/docs/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
% \input{"inst/docs/mqm/standard_example.txt"}
data(map10)                    # Mouse dataset

# simulate a cross (autosomes)
qtl <- c(3,15,3,7)             # QTL location
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01)

# MQM
crossaug <- mqmaugment(cross)  # Augmentation
cat(crossaug$mqm$Nind,'real individuals retained in dataset',
    crossaug$mqm$Naug,'individuals augmented')

result <- mqmscan(crossaug)    # Scan

# show LOD interval of the QTL on chr 3
lodint(result,chr=3)
% -----^^ inst/docs/mqm/standard_example.txt ^^-----
}

\keyword{utilities}
