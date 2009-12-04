\name{mqmaugment}
\alias{mqmaugment}

\title{MQM data augmentation}

\description{
The function $mqmaugment()$ fills in missing genotypes. For each missing marker
it fills in (augments) all possible genotypes and calculates their likelihood using
information on surrounding markers and recombination frequencies. When a
genotype is more likely than the $augment\_aboveprob$ parameter the
\emph{augmented} individual is stored in a new crossobject.
}

\usage{
mqmaugment(cross, pheno.col=1, maxaugind=60, augment_aboveprob=1, verbose=FALSE)
}

\arguments{
  \item{cross}{An object of class \code{cross}. See
   \code{\link[qtl]{read.cross}}.} 
  \item{pheno.col}{Column number in the phenotype matrix which should be
    used as the phenotype.  This can be a vector of integers, which index on
    the phenotype, a string matching one phenotype name, or, when the value is
    TRUE, all phenotypes are used.
  }
  \item{maxaugind}{ Maximum number of augmentations per individual }
  \item{augment_aboveprob}{ Augment genotypes that are above a probability of
  occurring (drop the genotypes with lower probabilities). Setting this value
  too high may result in dropping individuals entirely, thereby reducing the
  information in the dataset (using the VERBOSE switch will tell you).
  augment_aboveprob is a value between 0 and 1. For example a value of 0.5 will
  drop all genotypes that are half as likely than the most likely genotype. A
  value of 0.1 will drop all genotypes that are less likely of ocurring 1 in
  10, compared against the most likely genotype.
 }
  \item{verbose}{ Be verbose }
}

\details{
}

\value{
  The cross object the single phenotype \code{pheno.col}, with augmented
  individuals (many will be repeated multiple times). Individuals may have been
  dropped completely when all genotypes are below the augment_aboveprob. An
  added component to the cross object named \code{extra} contains information
  on which individuals are retained and repeated.  
}  

\author{R.C. Jansen ; modified by Pjotr Prins and Danny Arends}

\note{
The sex chromosome 'X' is dropped before augmentation.

Warning: with lots of missing data and/or chromosome markers a large
augment_aboveprob may cause the routine to drop a lot of indivuals from the
dataset. Run with verbose=TRUE to see how many individuals are retained in
the dataset.  

When augmentation goes out of bounds a warning is shown and the returned
should not be used.

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
}

\keyword{ utilities }
