\name{mqmscan}
\alias{mqmscan}

\title{Genome scan with a multiple QTL model (MQM)}

\description{
  Genome scan with a multiple QTL model.
}

\usage{
mqmscan(cross, cofactors=NULL, pheno.col = 1,
  model=c("additive","dominance"), forceML=FALSE,
  cofactor.significance=0.02, em.iter=1000,
  window.size=25.0, step.size=5.0,
  logtransform = FALSE, estimate.map = FALSE,
  plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
  multicore=TRUE, batchsize=10, n.clusters=1, test.normality=FALSE,off.end=0
  )
}

\arguments{
  \item{cross}{
An object of class \code{cross}. See \code{\link{read.cross}} for details. % \crossobject
  }
  \item{cofactors}{ List of cofactors to be analysed as cofactors in backward elimination
    procedure when building the QTL model. See \code{\link{mqmsetcofactors}} on how-to manually set cofactors
    for backward elimination. Or use \code{\link{mqmautocofactors}} for automatic selection of cofactors. Only
    three kind of (integer) values are allowed in the cofactor list. (0: no cofactor at this marker, 1: Use
    this marker as an additive cofactor, 2: Use this marker as an sexfactor (Dominant cofactor))
  }
  \item{pheno.col}{ Column number in the phenotype matrix which should be
    used as the phenotype.  This can be a vector of integers; One may also
    give a character strings matching the phenotype names.  Finally, one
    may give a numeric vector of phenotypeIDs. This should consist of
    integers with  0 < value < no. phenotypes.
  }
  \item{model}{ When scanning for QTLs should haplotype dominance be
    considered in an F2 intercross. Using the dominance model we scan
    for additive effects but also allow an additional effect where AA+AB
    versus BB and AA versus AB+BB. This setting is ignored for BC and RIL
    populations
  }
  \item{forceML}{ Specify which statistical method to use to estimate
    variance components to use when QTL modeling and mapping. Default
    usage is the Restricted maximum likelihood approach (REML). With this
    option a user can disable REML and use maximum likelihood.
  }
  \item{cofactor.significance}{ Significance level at which a cofactor is
    considered significant. This is estimated using an analysis of
    deviance, and compared to the level specified by the user. The
    cofactors that dont reach this level of statistical significance
    are NOT used in the mapping stage. Value between 0 and 1 }
  \item{em.iter}{ Maximum number of iterations for the EM algorithm to converge }
  \item{window.size}{ Window size for mapping QTL locations,
    this parameter is used in the interval mapping stage. When
    calculating LOD scores at a genomic position all cofactors
    within window.size are dropped to estimate the (unbiased)
    effect of the location under interest.
  }
  \item{step.size}{ Step size used in interval mapping. A lower
    step.size parameter increases the number of output points,
    this creates a smoother QTL profile
  }
  \item{off.end}{Distance (in cM) past the terminal markers on each
    chromosome to which the genotype simulations will be carried.
  }
  \item{logtransform}{ Indicate if the algorithm should do a log
    transformation on the trait data in the pheno.col
  }
  \item{estimate.map}{ Should Re-estimation of the marker locations
    on the genetic map occur before mapping QTLs. This method is
    deprecated rather use the \code{\link{est.map}} function in R/qtl.
    This is because no map is returned into the crossobject.
    The old map remains in the cross object.
    }
  \item{plot}{ plot the results (default FALSE)}
  \item{verbose}{ verbose output}
  \item{outputmarkers}{Needs to be explained}
  \item{multicore}{Use multicore (if available)}
  \item{batchsize}{Needs to be explained}
  \item{n.clusters}{Number of child processes to split the job into.}
  \item{test.normality}{If TRUE, test whether the phenotype follows a
    normal distribution via \code{\link{mqmtestnormal}}.
  }
}

%\details{
%}

\value{
When scanning a single phenotype the function returns a \code{\link{scanone}}
object.

The object contains a matrix of three columns for LOD scores, information content
and LOD*information content with pseudo markers sorted in increasing
order. For more information on the scanone object see: \code{\link{scanone}}
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\note{
The resulting scanone object itself can be visualized using the standard R/qtl
plotting routines (\code{\link{plot.scanone}}) or specialized function to show
the mqm model (\code{\link{mqmplot.singletrait}}) and QTL profile. If cofactors
were specified the QTL model used in scanning is also returned as a named
attribute of the scanone object called mqmmodel. It can be extracted from the
resulting scanone object by using the \code{\link{mqmgetmodel}} function or the
\code{\link{attr}} function.

Also note the \code{estimate.map} parameter does not return
its re-estimated genetic map, altough it is used internally.  When scanning
multiple genotypes a \code{mqmmulti} object is created.  This object is just a
list composed of scanone objects. The results for a single trait can be
obtained from the \code{mqmmulti} object, in scanone format.
}

\seealso{
  \itemize{
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

\keyword{ models }
