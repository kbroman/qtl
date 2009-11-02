\name{mqmaugment}
\alias{mqmaugment}

\title{MQM data augmentation}

\description{
Estimates missing data by increasing the number of marker locations and filling in 
genotypes based on surrounding markers and recombination frequencies. 

Augment/expand the dataset by adding additional (likely) genotypes.
Genotypic marker information is incomplete when marker scores are missing or dominant, or
markers are non-segregating but you want to assume a
segregating QTL on top of this marker (~very closely linked),
and all likely possible configurations are generated 

A new markerset is created, likewise the Augmentation halts when the number of
individuals maxNaug is reached. Individuals are expanded up to imaxNaug.
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
  occurring (ignore lower probabilities). Setting this value too high
  may result in dropping individuals entirely. E.g. for augment_aboveprob=100 : eliminate genotypes 100
  times less likely than the most likely configuration. The dataset may be
  reduced when there is a lot of missing data.  }
  \item{verbose}{ Be verbose }
}

\details{
}

\value{
  The cross object the single phenotype \code{pheno.col}, without dropped
  individuals and with augmented individuals (many will be repeated multiple
  times). An added component \code{extra} contains information on which
  individuals are retained and repeated.  
}  

\author{R.C. Jansen ; Danny Arends}

\note{
The sex chromosome 'X' is dropped before augmentation.

With unrealistic values for augmentation a warning is shown and the cross
object may have been corrupted. 

Warning: with lots of missing data and/or chromosome markers a neglect value
larger than 1 may cause the routine to drop a lot of indivuals from the
dataset. Try running with verbose to see how many individuals are retained in
the dataset.  
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
