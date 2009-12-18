\name{mqmscan}
\alias{mqmscan}
\alias{mqm}

\title{Genome scan with a multiple QTL model (MQM)}

\description{
  Genome scan with a multiple QTL model.
}

\usage{
mqm(cross, cofactors, pheno.col = 1, 
  model=c("Additive","Dominance"), forceML=false,
  cofactor.significance=0.02, em.iter=1000, 
  window.size=25.0, step.size=5.0, step.min=-20.0, step.max=max(unlist(pull.map(cross))),
  logtransform = FALSE, estimate.map = FALSE,
  plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
  multicore=TRUE, batchsize=10, n.clusters=1
  )
}

\arguments{
  \item{cross}{ 
    An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
  }
  \item{cofactors}{ List of cofactors to be analysed as cofactors in the 
    QTL model. See \code{\link{mqmsetcofactors}} on how-to set cofactors 
    for backward elimination. (0:no cofactor, 1:cofactor, 2:sexfactor) 
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
  \item{step.min}{ Begin this many cM before the chromosome, by default we 
    start at the genetic beginning (0cM) 
  }
  \item{step.max}{ Maximum mapping length in cM, by default set 
    to the maximum genetic length found in the crossobject. However
    when using the restimate.map parameter it could be needed
    to increase this parameter to accomodate for the larger map 
    internally estimated.
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
  \item{plot}{ plot the results (default TRUE)}
  \item{verbose}{ verbose output}
}

%\details{
%}

\value{
When scanning 1 phenotype the function returns a scanone object. This contains a matrix of three columns LODscores, 
Informationcontent and LOD*Informationcontent with at each row pseudomarkers
sorted in increasing order. For more information on the scanone object see: \code{\link{scanone}} 
The resulting scanone object itself can be visualized using the standard R/qtl 
plotting routines (\code{\link{plot.scanone}}) or specialized function to show
the mqm model (\code{\link{mqmplot_one}}) and QTL profile. If cofactors were 
specified the QTL model used in scanning is also returned as a named attribute 
of the scanone object called mqmmodel. It can be extracted from the resulting 
scanone object by using the \code{\link{mqmgetmodel}} function or the 
\code{\link{attr}} function. Finally the estimate.map paramter does not 
return its re estimated genetic map.
When scanning multiple genotypes a \code{mqmmulti} object is created.
This object is just a list composed of scanone objects. The results for a single
trait can be obtained from the \code{mqmmulti} object, in scanone format.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

%\note{ 
%}

\seealso{
\itemize{
% \input{"inst/docs/mqm/standard_seealso.txt"}
    \item \code{\link[qtl]{MQM}} - MQM description and references
    \item \code{\link[qtl]{mqmscan}} - Main MQM single trait analysis
    \item \code{\link[qtl]{mqmall}} - Parallellized traits analysis 
    \item \code{\link[qtl]{mqmaugment}} - Augmentation routine for estimating missing data
    \item \code{\link[qtl]{mqmcofactors}} - Set cofactors for the MQM algorithm
    \item \code{\link[qtl]{mqmbootstrap}} - Estimate confidence intervals
    \item \code{\link[qtl]{scanone}} - Single QTL scanning

% -----^^ inst/docs/mqm/standard_seealso.txt ^^-----
}
}

\examples{
% \input{"inst/docs/mqm/standard_example.txt"}
data(map10)                    # Mouse dataset
qtl <- c(3,15,3,7)             # QTL location

# simulate a cross (autosomes)
cross <- sim.cross(map10,qtl,n=100,missing.prob=0.01)

# MQM
crossaug <- mqmaugment(cross)  # Augmentation
result <- mqmscan(crossaug)    # Scan

# show LOD interval of the QTL on chr 3
lodint(result,chr=3)

cat(result$mqm$Nind,'real individuals retained in dataset', \
    result$mqm$Naug,'individuals augmented')
    result$mqm$Naug,'individuals augmented')

% -----^^ inst/docs/mqm/standard_example.txt ^^-----
}

\keyword{ models }
