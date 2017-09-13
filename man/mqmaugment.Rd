\name{mqmaugment}
\alias{mqmaugment}

\title{MQM augmentation}

\description{
Fill in missing genotypes for MQM mapping. For each missing or incomplete
marker it fills in (or `augments') all possible genotypes, thus creating new
candidate `individuals'. The probability of each indidual is calculated using
information on neighbouring markers and recombination frequencies. When a
genotype of an augmented genotype is less likely than the \code{minprob}
parameter it is dropped from the dataset. The \emph{augmented} list of
individuals is returned in a new cross object. For a full discussion on
augmentation see the MQM tutorial online.
}

\usage{
mqmaugment(cross, maxaugind=82, minprob=0.1,
           strategy=c("default","impute","drop"),
           verbose=FALSE)
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{maxaugind}{
    Maximum number of augmentations per individual.  The default of 82
    allows for six missing markers for an individual in a BC cross
    (\eqn{2^6=64}) and four missing markers in an F2 (\eqn{3^4=81}). When a
    large number of markers are missing this default number is quickly
    reached.
  }
  \item{minprob}{
    Return individuals with augmented genotypes that have at least this probability of
    occurring. \code{minprob} is a value between 0 and 1. For example a value of 0.5
    will drop all genotypes that are half as likely as the most likely
    genotype (candidate of the individual). The default value of 0.1 will drop
    all genotypes that are less likely of ocurring than 1 in 10, compared
    against the most likely genotype. Use a value of 1.0 to return a single
    filled in genotype for each individual.
 }
  \item{strategy}{
    When individuals have too much missing data and augmentation fails three
    options are provided:
    1. \code{"default"}: Calculate genotypes at missing marker positions,
    accounting for \code{minprob}, and add this individual to the set.
    2. \code{"impute"}: Calculate the most likely genotypes at missing marker
    positions and impute \code{maxaugind} individual-variants around the most
    likely genotype.
    3. \code{"drop"}: Drop individuals that cannot be augmented from the
    dataset, this option is not advised because information from the dropped
    individuals will be lost.
 }
  \item{verbose}{ If TRUE, give verbose output }
}

\value{
  Returns the cross object with augmented individuals (many individuals from
  the data set will be repeated multiple times). Some individuals may have been
  dropped completely when the probability falls below \code{minprob}. An added
  component to the cross object named \code{mqm} contains information on
  exactly which individuals are retained and repeated.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\note{
The sex chromosome 'X' is treated like autosomes during augmentation.
With an F2 the sex chromosome is not considered. This will change in
a future version of MQM.
Run with \code{verbose=TRUE} to verify how many individuals are augmented
versus moved to the second augmentation round. This could have an effect
on the resulting dataset or check the return \code{cross$mqm} values. Compare
results by using \code{minprob=1}.
}

\seealso{
  \itemize{
    \item \code{\link{fill.geno}} - Alternative routine for estimating missing data
% \input{"inst/doc/Sources/MQM/mqm/standard_seealso.txt"}
    \item The MQM tutorial: \url{https://rqtl.org/tutorials/MQM-tour.pdf}
    \item \code{\link{MQM}} - MQM description and references
    \item \code{\link{mqmscan}} - Main MQM single trait analysis
    \item \code{\link{mqmscanall}} - Parallellized traits analysis
    \item \code{\link{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link{mqmautocofactors}} - Set cofactors using marker density
    \item \code{\link{mqmsetcofactors}} - Set cofactors at fixed locations
    \item \code{\link{mqmpermutation}} - Estimate significance levels
    \item \code{\link{scanone}} - Single QTL scanning
% -----^^ inst/doc/Sources/MQM/mqm/standard_seealso.txt ^^-----
  }
}

\examples{
% \input{"inst/doc/Sources/MQM/mqm/standard_example.txt"}
data(map10)                    # Genetic map modeled after mouse

# simulate a cross (autosomes 1-10)
qtl <- c(3,15,1,0)             # QTL model: chr, pos'n, add've & dom effects
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01)

# MQM
crossaug <- mqmaugment(cross)  # Augmentation
cat(crossaug$mqm$Nind,'real individuals retained in dataset',
    crossaug$mqm$Naug,'individuals augmented\n')

result <- mqmscan(crossaug)    # Scan

# show LOD interval of the QTL on chr 3
lodint(result,chr=3)
% -----^^ inst/doc/Sources/MQM/mqm/standard_example.txt ^^-----
}

\keyword{utilities}
