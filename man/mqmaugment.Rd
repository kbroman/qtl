\name{mqmaugment}
\alias{mqmaugment}

\title{MQM data augmentation}

\description{
Fills in missing genotypes for multiple QTL model (MQM) mapping. For each
missing marker it fills in (or 'augments') all possible genotypes, thus
creating new candidate 'individuals' and calculates their probability using
information on neighbouring markers and recombination frequencies. When a
genotype of an augmented genotype is less likely than the $minprob$ parameter
it is dropped from the dataset. The \emph{augmented} list of individuals is
returned in a new $crossobject$.  
}

\usage{
mqmaugment(cross, pheno.col=1, maxaugind=60, minprob=0.1, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}}.} 
  \item{pheno.col}{Column number in the phenotype matrix which should be
    used as the phenotype. This can be a vector of integers, which index on
    the phenotype, a string matching one phenotype name, or, when the value is
    TRUE, all phenotypes are used. (Note: the current implementation
    allows one phenotype only).
  }
  \item{maxaugind}{ Maximum number of augmentations per individual.
  The default of 60 allows for 20 missing markers for an individual in
  a BC cross and 30 missing markers in an F2. When a large number of
  markers are missing this default number is quickly reached.
  }
  \item{minprob}{ 
    Return individuals with augmented genotypes that are above a
    probability of occurring.  Setting this value too low may result
    in dropping individuals entirely, thereby reducing the information
    in the dataset (using the VERBOSE switch will tell you). minprob
    is a value between 0 and 1. For example a value of 0.5 will drop
    all genotypes that are half as likely than the most likely
    genotype (individual). The default value of 0.1 will drop all
    genotypes that are less likely of ocurring than 1 in 10, compared
    against the most likely genotype. Use a value of 1 to fill in only the
	most likely genotype (no augmentation, no dropping of individuals). \FIXME
 }
  \item{verbose}{ Verbose output }
}

\details{
}

\value{
  Returns the cross object for the single phenotype \code{pheno.col},
  with augmented individuals (many individuals from the data set will
  be repeated multiple times). Some individuals may have been dropped
  completely when the probability falls below $minprob$. \FIXME An added
  component to the cross object named \code{extra} contains
  information on exactly which individuals are retained and repeated.  
}  

\author{Original algorithm by R.C. Jansen ; modified by Danny Arends
and Pjotr Prins}

\note{
The sex chromosome 'X' is treated like autosomes during augmentation.
With an F2 the sex chromosome is not considered. This will change in
a future version of MQM.

When augmentation runs out of memory, a warning is given and the returned
object should not be used.

With a large amount of missing data a high\FIXME $minprob$ may cause the
routine to drop a lot of indivuals from the dataset. Run with
verbose=TRUE to verify how many individuals are retained in the
dataset or check the return cross$extra values.
}

\seealso{
  \itemize{
  \item \code{\link[qtl]{MQM}} - MQM description and references
  \item \code{\link[qtl]{mqmscan}} - Main MQM function for trait analysis
  }
}

\examples{
qtl <- c(3,15,3,7)                     # QTL at chromosome 3
data(map10)                            # Mouse genome
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01)    # Simulate a Cross (just autosomes)
crossaug <- mqmaugment(cross)          # Augmentation of the data
result <- mqmscan(crossaug)            # Do a MQM scan of the genome
cat(result$extra$Nind,'retained individuals',result$extra$Naug,'augmented individuals')

}

\keyword{utilities}
