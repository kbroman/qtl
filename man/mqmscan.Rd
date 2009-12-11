\name{mqmscan}
\alias{mqmscan}
\alias{mqm}

\title{Genome scan with a multiple QTL model (MQM)}

\description{
  Genome scan with a multiple QTL model.
}

\usage{
	mqmscan(cross, cofactors, pheno.col = 1, REMLorML = 0, 
          alfa = 0.02, em.iter = 1000, windowsize = 25.0, step.size = 5.0, 
          step.min = -20.0, step.max = 220.0, file = "MQM_output.txt",
          doLOG=0, est.map=0, dominance=0, plot=FALSE, forceRIL=0, 
          verbose=FALSE)
}

\arguments{
  \item{cross}{ 
An object of class \code{cross}. See \code{\link{read.cross}} for details % \crossobject
    }
  \item{cofactors}{ List of cofactors with length marker (0:no cofactor, 1:cofactor, 2:sexfactor) }
  \item{pheno.col}{ Column number of the phenotype under analysis }    
  \item{REMLorML}{ Use reduced maximum likelihood (0) or maximum likelihood (1) }
  \item{alfa}{ Global p-value (0<= alfa <= 1) }
  \item{em.iter}{ Maximum number of iterations EM algorithm }
  \item{windowsize}{ Window size for mapping QTL locations }
  \item{step.size}{ Step size used for the generation of output QTL }
  \item{step.min}{ Begin this many cM before the chromosome }
  \item{step.max}{ Maximum length in cM }
  \item{file}{ File for raw results}
  \item{doLOG}{ Indicate if the algorithm should LOG transform traitdata (DEFAULT 0 = off | 1 = on) }    
  \item{est.map}{ Reestimation of the marker locations on the genetic map before mapping QTLs (no map is returned the old map remains in the cross object }
  \item{dominance}{ Should dominance be taken into account ( 0= Additive  (DEFAULT)/ 1= Dominance) }  
  \item{forceRIL}{ Should we force the cross as a RIL, this option removed AB genotypes and replaces them with BB ? ( 0= NO  (DEFAULT)/ 1= YES) }  
  \item{plot}{ plot the results (default TRUE)}
  \item{verbose}{ verbose output}
}

\details{
}

\value{
Returns a matrix with QTLscores,Informationcontent,QTL*Informationcontent for
each location on the chromosome, which can be plotted using the plot.scanone
function of R/QTL or plot.MQMone to also visualize the QTL*informationcontent
The phenotype under consideration is specified by the pheno.col argument of the
function.
}

\author{
Ritsert C Jansen; Danny Arends; Pjotr Prins; Karl W Broman \email{kbroman@biostat.wisc.edu} % \mqmauthors
}

\note{ 
}

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
cross <- sim.cross(map10[1:10],qtl,n=100,missing.prob=0.01)

# MQM
crossaug <- mqmaugment(cross)  # Augmentation
result <- mqmscan(crossaug)    # Scan

# show LOD interval of the QTL on chr 3
lodint(result,chr=3)

cat(result$extra$Nind,'real individuals retained in dataset', \
    result$extra$Naug,'individuals augmented')

% -----^^ inst/docs/mqm/standard_example.txt ^^-----
}

\keyword{ models }
